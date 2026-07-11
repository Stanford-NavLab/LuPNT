"""Authoring a new Application entirely in Python.

Demonstrates pylupnt's Python-authoring hooks: subclass ``pnt.Application``, register it with
``pnt.register_application(name, cls)``, and reference it by ``class:`` in a scenario config --
the C++ ``pnt.Simulation`` instantiates and drives it. Here the app is a full angles-only
orbit-determination EKF (numpy), whose filter dynamics are a pylupnt ``NBodyDynamics`` built and
propagated (with STM) from Python via ``dyn.propagate_stm``. The observer estimates its own orbit
from unit line-of-sight bearings to a known-ephemeris target; results are stored on ``self`` and
read back after ``sim.run()``. Converges to a few metres over the 6 h arc.

Run:  pixi run python python/examples/py_authored_app_demo.py
"""

import sys, os
import numpy as np
sys.path.insert(0, "python"); os.environ["LUPNT_DATA_PATH"] = os.path.abspath("data/LuPNT_data")
import yaml, pylupnt as pnt


def make_dyn():
    d = pnt.NBodyDynamics()
    d.add_body(pnt.Body.Moon(8, 8)); d.add_body(pnt.Body.Earth()); d.add_body(pnt.Body.Sun())
    d.set_frame(pnt.Frame.MOON_CI); d.set_integrator(pnt.IntegratorType.RKF45); d.set_autodiff(True)
    return d


class PyAnglesOdtsApp(pnt.Application):
    def __init__(self, config):
        pnt.Application.__init__(self)
        self.target = config["target"]
        self.set_frequency(config.get("frequency", 1.0 / 60))
        self.sigma = np.radians(config.get("angle_sigma_arcsec", 2.0) / 3600.0)
        self.p0 = config.get("initial_position_sigma_m", 300.0)
        self.v0 = config.get("initial_velocity_sigma_mps", 0.3)
        self.qa = config.get("process_accel_sigma_mps2", 1e-6)
        self.rng = np.random.default_rng(config.get("seed", 42))
        self.times, self.est, self.truth, self.sig = [], [], [], []
        self._init = False

    def setup(self):
        pnt.Application.setup(self)
        self.dyn = make_dyn()

    def step(self, t):
        if not self._init:                                   # capture initial truths once
            self.obs = np.asarray(self.get_agent().get_state_at(0.0))
            self.tgt = np.asarray(self.get_agent().get_world().get_state_at(self.target, 0.0))
            self.x = self.obs + np.r_[self.rng.normal(0, self.p0, 3), self.rng.normal(0, self.v0, 3)]
            self.P = np.diag(np.r_[[self.p0**2]*3, [self.v0**2]*3])
            self._init, self._t = True, t
            self._log(t); return
        dt = t - self._t
        # propagate the truth references and the estimate consistently with the same model
        self.obs = np.asarray(self.dyn.propagate_stm(self.obs, float(self._t), float(t))[0]).ravel()
        self.tgt = np.asarray(self.dyn.propagate_stm(self.tgt, float(self._t), float(t))[0]).ravel()
        xf, F = self.dyn.propagate_stm(self.x, float(self._t), float(t))
        self.x = np.asarray(xf).ravel(); F = np.asarray(F); self._t = t
        q = self.qa**2
        Q = q * np.block([[dt**3 / 3 * np.eye(3), dt**2 / 2 * np.eye(3)],
                          [dt**2 / 2 * np.eye(3), dt * np.eye(3)]])
        self.P = F @ self.P @ F.T + Q
        r = self.x[:3]; d = self.tgt[:3] - r; rn = np.linalg.norm(d); u = d / rn
        H = np.zeros((3, 6)); H[:, :3] = -(np.eye(3) - np.outer(u, u)) / rn
        u_t = self.tgt[:3] - self.obs[:3]; u_t = u_t / np.linalg.norm(u_t)
        z = u_t + self.rng.normal(0, self.sigma, 3)
        R = self.sigma**2 * np.eye(3)
        S = H @ self.P @ H.T + R
        K = self.P @ H.T @ np.linalg.inv(S)
        self.x = self.x + K @ (z - u)
        IKH = np.eye(6) - K @ H
        self.P = IKH @ self.P @ IKH.T + K @ R @ K.T
        self._log(t)

    def _log(self, t):
        self.times.append(t); self.est.append(self.x.copy())
        self.truth.append(self.obs.copy()); self.sig.append(np.sqrt(np.diag(self.P)))


pnt.register_application("PyAnglesOdtsApp", PyAnglesOdtsApp)
cfg = yaml.safe_load(open("configs/sat_bearing_odts.yaml"))
a = dict(cfg["agents"]["observer"]["application"]); a["class"] = "PyAnglesOdtsApp"
a["process_accel_sigma_mps2"] = 1e-6; a.pop("monte_carlo_runs", None)
cfg["agents"]["observer"]["application"] = a
sim = pnt.Simulation(cfg); sim.run()
app = sim.get_agent("observer").get_application()
est = np.array(app.est); truth = np.array(app.truth)
perr = np.linalg.norm(est[:, :3] - truth[:, :3], axis=1)
print(f"PY-EKF steps={len(perr)} init={perr[0]:.1f} final={perr[-1]:.1f} "
      f"RMS(last20%)={np.sqrt(np.mean(perr[int(0.8*len(perr)):]**2)):.1f} m")
