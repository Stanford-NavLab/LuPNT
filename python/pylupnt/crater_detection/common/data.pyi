from __future__ import annotations
import datetime as dt
import h5py as h5py
from matplotlib import pyplot as plt
import numpy as np
import os as os
import pylupnt.crater_detection.common.conics
from pylupnt.crater_detection.common.conics import MaskGenerator
from pylupnt.crater_detection.common import constants as const
from tqdm.asyncio import tqdm_asyncio as tq
import uuid as uuid
__all__ = ['DataGenerator', 'MaskGenerator', 'const', 'demo_settings', 'dt', 'generate', 'h5py', 'inspect_dataset', 'make_dataset', 'np', 'os', 'plt', 'tq', 'uuid']
class DataGenerator(pylupnt.crater_detection.common.conics.MaskGenerator):
    def __init__(self, *args, **kwargs):
        ...
    def image_mask_pair(self, **mask_kwargs):
        ...
def demo_settings(n_demo = 20, generation_kwargs = None):
    ...
def generate(size, **kwargs):
    ...
def inspect_dataset(dataset_path, plot = True, summary = True, n_inspect = 25, pixel_range = (0, 1), return_fig = False):
    ...
def make_dataset(n_training, n_validation, n_testing, output_path = None, identifier = None, generation_kwargs = None):
    ...
