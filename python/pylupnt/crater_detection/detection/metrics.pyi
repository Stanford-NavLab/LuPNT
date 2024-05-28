from __future__ import annotations
import numpy as np
from numpy import linalg as LA
from pylupnt.crater_detection.common.conics import conic_center
from pylupnt.crater_detection.common.conics import ellipse_axes
from pylupnt.crater_detection.common.conics import scale_det
import torch as torch
from torchvision.ops.boxes import box_iou
__all__ = ['LA', 'box_iou', 'conic_center', 'detection_metrics', 'ellipse_axes', 'f1_score', 'gaussian_angle_distance', 'get_matched_idxs', 'mv_kullback_leibler_divergence', 'norm_mv_kullback_leibler_divergence', 'np', 'precision_recall', 'scale_det', 'torch']
def detection_metrics(pred_dict: typing.Dict, target_dict: typing.Dict, iou_threshold: float = 0.5, confidence_threshold: float = 0.75, distance_threshold: float = None) -> typing.Tuple[float, float, float, float, float]:
    """
    
        Calculates Precision, Recall, F1, IoU, and Gaussian Angle Distance for a single image.
    
        Parameters
        ----------
        pred_dict
            Prediction output dictionary (from CraterDetector)
        target_dict
            Target dictionary
        iou_threshold
            Minimum IoU to consider prediction and target to be matched
        confidence_threshold
            Minimum prediction class score
        distance_threshold
            (default: None) If given, use to add another condition for determining TP, FP, FN
    
        Returns
        -------
        Precision, Recall, F1, IoU, Gaussian Angle Distance
        
    """
def f1_score(precision: float, recall: float) -> float:
    """
    
        Calculates F1 score according to:
    
        .. math:: 2 (P*R)/(P+R)
    
        Parameters
        ----------
        precision
            Precision
        recall
            Recall
    
        Returns
        -------
        F1 score
    
        
    """
def gaussian_angle_distance(A1: typing.Union[torch.Tensor, numpy.ndarray], A2: typing.Union[torch.Tensor, numpy.ndarray]) -> torch.Tensor:
    ...
def get_matched_idxs(pred: typing.Union[typing.Dict, torch.Tensor], target: typing.Union[typing.Dict, torch.Tensor], iou_threshold: float = 0.5, return_iou: bool = False) -> typing.Tuple:
    """
    
        Returns indices at which IoU is maximum, as well as a mask containing whether it's above iou_threshold.
    
        Parameters
        ----------
        pred
            Prediction boxes or dictionary
        target
            Target bounding boxes or dictionary
        iou_threshold
            Minimum IoU to consider a prediction a True Positive
        return_iou
            Whether to return IoU values of matched detections
    
        Returns
        -------
        Matching indices, boolean match mask
    
        Examples
        --------
        >>> matched_idxs, matched, iou_list = get_matched_idxs(boxes_pred, boxes_target, return_iou=True)
        >>> pred_true = pred[matched]
        >>> target_matched = target[matched_idxs][matched]
    
        
    """
def mv_kullback_leibler_divergence(A1: torch.Tensor, A2: torch.Tensor, shape_only: bool = False) -> torch.Tensor:
    ...
def norm_mv_kullback_leibler_divergence(A1: torch.Tensor, A2: torch.Tensor) -> torch.Tensor:
    ...
def precision_recall(TP: int, FP: int, FN: int) -> typing.Tuple[float, float]:
    """
    
        Calculates Precision and Recall from detections according to:
    
        .. math:: P = TP / (TP + FP)
        .. math:: R = TP / (TP + FN)
    
        Parameters
        ----------
        TP
            True Positives
        FP
            False Positives
        FN
            False Negatives
    
        Returns
        -------
        Precision, Recall
    
        
    """
