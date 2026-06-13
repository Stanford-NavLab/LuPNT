# Plasmasphere Delay Data Generation

This folder provides a set of python scripts and jupyter notebooks to generate GNSS measurement data with simulated plasmasphere delays

## Pre-requisits

## Data Generation Steps
Recommended workflow is as below. The process is divided into multiple steps since each process requires a lot of compute time.

1. Run `generate_labels.py` to generate pickle files of gnss data
2. Run `convert_labels_to_csv.py` to convert the pickle files into csv
3. Run `run_raytrace_batch.py` with `freq_families = [1]` and `straight_rays = [True]` to compute the raytracing in c++ (for L1/E1 signals)
4. Run `convert_L1_to_L5.py` to convert the delays for L1 to L5
5. Run `run_raytrace.py` again with `freq_families = [5]` and `straight_rays = [True]` to simulate L5/E5a signals that are not covered by convergion
6. (Optional) The ray-tracing results for 3 and 5 does not account for ray-bending. This is fine for most of the signals ( with tangential above 2000 km Earth), but for higher fidelity simulation, it will be nice to simulate ray-bending for lower altitude rays. For this, rerun `run_raytrace.py` with `straight_rays = [False]`
7. Run `convert_ionodata.ipynb` to convert the raytracing results to files to pickle and csv files.
8. Run `generate_orbit_h5.ipynb` to convert the propagated orbits into h5 format.

## Simulation Parameters
The major parameters to change for data generation is follows
- Starting Epoch: change `gps_datetime` in `generate_labels.py`
- Simulation Time Lenght: change `n_orbit ` in `generate_labels.py`
- Measurement time inverval: change `dt` in `generate_labels.py`
- User index: The receiver (user) can be either of the 5 LCRNS satellites or fixed user on lunar south-pole. Change `lcrns_idxs` in `generate_labels.py`
- Rz12:  The index for solar activity, change `rz12` in `run_raytrace.py`
- kp: The index for geomagnetic activity, change `kp` in `run_raytrace.py`
