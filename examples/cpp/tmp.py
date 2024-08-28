import pylupnt as pnt
import numpy as np

body = pnt.Body.Moon(50, 50)

# Dynamics
dyn_nbody50 = pnt.NBodyDynamics(pnt.IntegratorType.RK4)
dyn_nbody50.add_body(pnt.Body.Moon(50, 50))
dyn_nbody50.add_body(pnt.Body.Earth())
dyn_nbody50.add_body(pnt.Body.Sun())
dyn_nbody50.set_time_step(60)
dyn_nbody50.set_frame(pnt.MOON_CI)


t0 = 300891599.99999976
tfs = np.array([300891599.99999976, 300895199.99999976, 300898799.99999976])
rv0_mi = np.array(
    [
        1076.8119928780047,
        -334.33334391607241,
        2361.1614474109979,
        -0.87564912639352988,
        1.3707947272813383,
        0.59344178577225781,
    ]
)
rv_case5_mi = dyn_nbody50.propagate(rv0_mi, t0, tfs, progress=True)
print(rv_case5_mi)
