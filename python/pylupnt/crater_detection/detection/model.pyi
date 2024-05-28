from __future__ import annotations
from pylupnt.crater_detection.detection.roi_heads import EllipseRegressor
from pylupnt.crater_detection.detection.roi_heads import EllipseRoIHeads
import torch as torch
from torch.nn.modules.conv import Conv2d
from torchvision.models.detection.anchor_utils import AnchorGenerator
from torchvision.models.detection.backbone_utils import resnet_fpn_backbone
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor
from torchvision.models.detection.faster_rcnn import TwoMLPHead
import torchvision.models.detection.generalized_rcnn
from torchvision.models.detection.generalized_rcnn import GeneralizedRCNN
from torchvision.models.detection.rpn import RPNHead
from torchvision.models.detection.rpn import RegionProposalNetwork
from torchvision.models.detection.transform import GeneralizedRCNNTransform
from torchvision.ops.poolers import MultiScaleRoIAlign
__all__ = ['AnchorGenerator', 'Conv2d', 'CraterDetector', 'EllipseRegressor', 'EllipseRoIHeads', 'FastRCNNPredictor', 'GeneralizedRCNN', 'GeneralizedRCNNTransform', 'MultiScaleRoIAlign', 'RPNHead', 'RegionProposalNetwork', 'TwoMLPHead', 'resnet_fpn_backbone', 'torch']
class CraterDetector(torchvision.models.detection.generalized_rcnn.GeneralizedRCNN):
    @staticmethod
    def get_conics(*args, **kwargs):
        ...
    def __init__(self, num_classes = 2, backbone_name = 'resnet50', min_size = 256, max_size = 512, image_mean = None, image_std = None, rpn_anchor_generator = None, rpn_head = None, rpn_pre_nms_top_n_train = 2000, rpn_pre_nms_top_n_test = 1000, rpn_post_nms_top_n_train = 2000, rpn_post_nms_top_n_test = 1000, rpn_nms_thresh = 0.7, rpn_fg_iou_thresh = 0.7, rpn_bg_iou_thresh = 0.3, rpn_batch_size_per_image = 256, rpn_positive_fraction = 0.5, rpn_score_thresh = 0.0, box_roi_pool = None, box_head = None, box_predictor = None, box_score_thresh = 0.05, box_nms_thresh = 0.5, box_detections_per_img = 100, box_fg_iou_thresh = 0.5, box_bg_iou_thresh = 0.5, box_batch_size_per_image = 512, box_positive_fraction = 0.25, bbox_reg_weights = None, ellipse_roi_pool = None, ellipse_head = None, ellipse_predictor = None, ellipse_loss_metric = 'gaussian-angle'):
        ...
