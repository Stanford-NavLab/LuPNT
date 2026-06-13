import torch
import torch.nn as nn
import torch.nn.functional as F
from core.update import BasicMultiUpdateBlock, LSTMMultiUpdateBlock
from core.extractor import BasicEncoder, MultiBasicEncoder, ResidualBlock
from core.corr import CorrBlock1D, PytorchAlternateCorrBlock1D, CorrBlockFast1D, AlternateCorrBlock
from core.utils.utils import coords_grid, upflow
from nets.refinement import *

try:
    autocast = torch.cuda.amp.autocast
except:
    # dummy autocast for PyTorch < 1.6
    class autocast:
        def __init__(self, enabled):
            pass

        def __enter__(self):
            pass

        def __exit__(self, *args):
            pass


class EAIStereo(nn.Module):
    def __init__(self, args):
        super().__init__()
        self.args = args

        context_dims = args.hidden_dims

        print('args.hidden_dims')
        print(len(args.hidden_dims))

        self.cnet = MultiBasicEncoder(output_dim=[args.hidden_dims, context_dims], norm_fn="batch",
                                      downsample=args.n_downsample)  #
        self.update_block = LSTMMultiUpdateBlock(self.args, hidden_dims=args.hidden_dims)

        self.mh_convs = nn.ModuleList(
            [nn.Conv2d(context_dims[i], args.hidden_dims[i] * 4, 3, padding=3 // 2) for i in
             range(self.args.n_gru_layers)])

        self.conv2 = nn.Sequential(
            ResidualBlock(128, 128, 'instance', stride=1),
            nn.Conv2d(128, 256, 3, padding=1))

        self.errorAwareRefinement = errorAwareRefinement()

    def freeze_bn(self):
        for m in self.modules():
            if isinstance(m, nn.BatchNorm2d):
                m.eval()

    def initialize_flow(self, img):
        """ Flow is represented as difference between two coordinate grids flow = coords1 - coords0"""
        N, _, H, W = img.shape

        coords0 = coords_grid(N, H, W).to(img.device)
        coords1 = coords_grid(N, H, W).to(img.device)

        return coords0, coords1

    def upsample_flow(self, flow, mask):
        """ Upsample flow field [H/8, W/8, 2] -> [H, W, 2] using convex combination """
        N, D, H, W = flow.shape
        factor = 2 ** self.args.n_downsample
        mask = mask.view(N, 1, 9, factor, factor, H, W)
        mask = torch.softmax(mask, dim=2)

        up_flow = F.unfold(factor * flow, [3, 3], padding=1)
        up_flow = up_flow.view(N, D, 9, 1, 1, H, W)

        up_flow = torch.sum(mask * up_flow, dim=2)
        up_flow = up_flow.permute(0, 1, 4, 2, 5, 3)
        return up_flow.reshape(N, D, factor * H, factor * W)

    def forward(self, image1, image2, iters=12, flow_init=None, test_mode=False):
        """ Estimate optical flow between pair of frames """

        image1 = (2 * (image1 / 255.0) - 1.0).contiguous()
        image2 = (2 * (image2 / 255.0) - 1.0).contiguous()

        # run the context network
        with autocast(enabled=self.args.mixed_precision):
            *cnet_list, x = self.cnet(torch.cat((image1, image2), dim=0), dual_inp=True,
                                      num_layers=self.args.n_gru_layers)
            fmap1, fmap2 = self.conv2(x).split(dim=0, split_size=x.shape[0] // 2)
            net_h = [torch.tanh(x[0]) for x in cnet_list]  #
            net_ext = [torch.relu(x[1]) for x in cnet_list]  #

            net_ext = [list(conv(i).split(split_size=conv.out_channels // 4, dim=1)) for i, conv in
                       zip(net_ext, self.mh_convs)]


        if self.args.corr_implementation == "reg":  # Default
            corr_block = CorrBlock1D
            fmap1, fmap2 = fmap1.float(), fmap2.float()
        elif self.args.corr_implementation == "reg_cuda":  # Faster version of reg
            corr_block = CorrBlockFast1D
        corr_fn = corr_block(fmap1, fmap2, radius=self.args.corr_radius, num_levels=self.args.corr_levels)

        coords0, coords1 = self.initialize_flow(net_h[0])  # Context Encoder ������ڳ�ʼ���Ӳ�ͼ

        if flow_init is not None:
            coords1 = coords1 + flow_init

        flow_predictions = []
        cnt = 0
        for itr in range(iters):
            coords1 = coords1.detach()
            corr = corr_fn(coords1)  #
            flow = coords1 - coords0
            with autocast(enabled=self.args.mixed_precision):
                if self.args.n_gru_layers == 3 and self.args.slow_fast_gru:  # Update low-res GRU
                    net_h = self.update_block(net_h, net_ext, iter32=True, iter16=False, iter08=False,
                                              update=False)
                if self.args.n_gru_layers >= 2 and self.args.slow_fast_gru:  # Update low-res GRU and mid-res GRU
                    net_h = self.update_block(net_h, net_ext, iter32=self.args.n_gru_layers == 3, iter16=True,
                                              iter08=False, update=False)
                if cnt == 0:
                    netC = net_h
                    cnt = 1
                netC, net_h, up_mask, delta_flow = self.update_block(netC, net_h, net_ext, corr, flow,
                                                                     iter32=self.args.n_gru_layers == 3,
                                                                     iter16=self.args.n_gru_layers >= 2)

            delta_flow[:, 1] = 0.0

            coords1 = coords1 + delta_flow

            if test_mode and itr < iters - 1:
                continue

            if up_mask is None:
                disp_fullres = upflow(coords1 - coords0)
            else:
                disp_fullres = self.upsample_flow(coords1 - coords0, up_mask)
            disp_fullres = disp_fullres[:, :1]

            if itr == iters - 1:
                # disp_fullres[disp_fullres < 0] = 0
                if disp_fullres.max() < 0:
                    refine_value = self.errorAwareRefinement(disp_fullres, image1, image2)
                    disp_fullres = disp_fullres + refine_value
                else:
                    pass
                    # print('initializing...')

            # print(disp_fullres)

            flow_predictions.append(disp_fullres)

        if test_mode:
            return coords1 - coords0, disp_fullres

        return flow_predictions
