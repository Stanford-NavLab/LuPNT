# EAI-Stereo

## Software Requirements
PyTorch 1.12.0\
CUDA 11.7

```Shell
pip install scipy
pip install tqdm
pip install tensorboard
pip install opt_einsum
pip install imageio
pip install opencv-python
pip install scikit-image
```

## Required Data
To evaluate/train EAI-Stereo, you will need to download the required datasets. 
* [Sceneflow](https://lmb.informatik.uni-freiburg.de/resources/datasets/SceneFlowDatasets.en.html#:~:text=on%20Academic%20Torrents-,FlyingThings3D,-Driving) 
* [Middlebury](https://vision.middlebury.edu/stereo/data/)
* [ETH3D](https://www.eth3d.net/datasets#low-res-two-view-test-data)
* [KITTI](http://www.cvlibs.net/datasets/kitti/eval_scene_flow.php?benchmark=stereo)

By default `stereo_datasets.py` will search for the datasets in these locations. You can create symbolic links to wherever the datasets were downloaded in the `datasets` folder

```Shell
├── datasets
    ├── FlyingThings3D
        ├── frames_cleanpass
        ├── frames_finalpass
        ├── disparity
    ├── Monkaa
        ├── frames_cleanpass
        ├── frames_finalpass
        ├── disparity
    ├── Driving
        ├── frames_cleanpass
        ├── frames_finalpass
        ├── disparity
    ├── KITTI
        ├── testing
        ├── training
        ├── devkit
    ├── Middlebury
        ├── MiddEval3
    ├── ETH3D
        ├── two_view_testing
```


## Build sampler
```Shell
cd sampler && python setup.py install && cd ..
```

## Train
```Shell
bash ./train.sh
```

## Evaluate
Set the arguments in evaluate_stereo.py and execute
```Shell
python evaluate_stereo.py
```
