from __future__ import annotations
import numpy as np
from pylupnt.crater_detection.common.conics import conic_matrix
from pylupnt.crater_detection.detection.metrics import gaussian_angle_distance
from pylupnt.crater_detection.detection.metrics import norm_mv_kullback_leibler_divergence
import torch as torch
from torch import nn
import torchvision.models.detection.roi_heads
from torchvision.models.detection.roi_heads import RoIHeads
from torchvision.models.detection.roi_heads import fastrcnn_loss
__all__ = ['EllipseRegressor', 'EllipseRoIHeads', 'RoIHeads', 'conic_matrix', 'ellipse_loss_GA', 'ellipse_loss_KLD', 'fastrcnn_loss', 'gaussian_angle_distance', 'nn', 'norm_mv_kullback_leibler_divergence', 'np', 'postprocess_ellipse_predictor', 'preprocess_ellipse_targets', 'torch']
class EllipseRegressor(torch.nn.modules.module.Module):
    def __init__(self, in_channels = 1024, hidden_size = 512, out_features = 3):
        ...
    def forward(self, x) -> torch.Tensor:
        ...
class EllipseRoIHeads(torchvision.models.detection.roi_heads.RoIHeads):
    def __init__(self, box_roi_pool, box_head, box_predictor, fg_iou_thresh, bg_iou_thresh, batch_size_per_image, positive_fraction, bbox_reg_weights, score_thresh, nms_thresh, detections_per_img, ellipse_roi_pool, ellipse_head, ellipse_predictor, ellipse_loss_metric = 'gaussian-angle'):
        ...
    def forward(self, features: typing.Dict[str, torch.Tensor], proposals: typing.List[torch.Tensor], image_shapes: typing.List[typing.Tuple[int, int]], targets: typing.Optional[typing.List[typing.Dict[str, torch.Tensor]]] = None) -> typing.Tuple[typing.List[typing.Dict[str, torch.Tensor]], typing.Dict[str, torch.Tensor]]:
        ...
    def has_ellipse_reg(self) -> bool:
        ...
def ellipse_loss_GA(d_pred: torch.Tensor, ellipse_matrix_targets: typing.List[torch.Tensor], pos_matched_idxs: typing.List[torch.Tensor], boxes: typing.List[torch.Tensor]) -> torch.Tensor:
    ...
def ellipse_loss_KLD(d_pred: torch.Tensor, ellipse_matrix_targets: typing.List[torch.Tensor], pos_matched_idxs: typing.List[torch.Tensor], boxes: typing.List[torch.Tensor]) -> torch.Tensor:
    ...
def postprocess_ellipse_predictor(d_a: torch.Tensor, d_b: torch.Tensor, d_angle: torch.Tensor, boxes: torch.Tensor) -> torch.Tensor:
    ...
def preprocess_ellipse_targets(A_target: torch.Tensor):
    ...
