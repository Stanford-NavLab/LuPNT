import numpy as np


class NavSatFailureModel:
    """
    Piecewise 'bathtub' model:
      [0, tau1): Weibull (beta0, eta0) with target S(tau1)=S_tau1
      [tau1, tau2): Exponential with rate lam1
      [tau2, inf): Weibull (beta2, eta2) fitted to hit S targets at t12, t15
    """

    def __init__(
        self,
        tau1=0.25,  # years (≈3 months)
        S_tau1=0.985,  # survival at tau1
        beta0=0.6,  # infant mortality shape
        lam1=0.008,  # /year, nominal constant hazard
        tau2=10.0,  # start of wear-out (years)
        targets_years=(12.0, 15.0),  # targets for wear-out survival
        targets_survival=(0.88, 0.75),  # S(t12), S(t15)
    ):
        self.tau1, self.S_tau1, self.beta0 = tau1, S_tau1, beta0
        self.lam1, self.tau2 = lam1, tau2

        # Solve infant-mortality scale eta0 from S(tau1)
        # S(t) = exp(-(t/eta0)^beta0) => eta0 = tau1 / (-ln S_tau1)^(1/beta0)
        self.eta0 = tau1 / (-np.log(S_tau1)) ** (1.0 / beta0)

        # Survival at tau2 (end of nominal phase)
        self.S_tau2 = S_tau1 * np.exp(-lam1 * (tau2 - tau1))

        # Fit wear-out Weibull (beta2, eta2) to hit S(12)=S12 and S(15)=S15
        # For t>=tau2: S(t) = S_tau2 * exp(-((t-tau2)/eta2)^beta2)
        # Let A12 = -ln(S12/S_tau2) = ( (12-tau2)/eta2 )^beta2
        # and A15 = -ln(S15/S_tau2) = ( (15-tau2)/eta2 )^beta2
        # => (A12/A15) = ((12-tau2)/(15-tau2))^beta2  -> solve beta2
        dt1 = targets_years[0] - tau2
        dt2 = targets_years[1] - tau2

        A1 = -np.log(targets_survival[0] / self.S_tau2)
        A2 = -np.log(targets_survival[1] / self.S_tau2)
        # guard for numeric sanity
        if not (A1 > 0 and A2 > 0 and dt1 > 0 and dt2 > 0):
            raise ValueError("Targets produce invalid wear-out fit.")
        self.beta2 = np.log(A1 / A2) / np.log(dt1 / dt2)
        self.eta2 = dt1 / (A1 ** (1.0 / self.beta2))

    def survival(self, t):
        t = np.asarray(t, dtype=float)
        S = np.ones_like(t)
        # Region 1: 0 <= t < tau1
        m1 = (t >= 0) & (t < self.tau1)
        S[m1] = np.exp(-((t[m1] / self.eta0) ** self.beta0))

        # Region 2: tau1 <= t < tau2
        m2 = (t >= self.tau1) & (t < self.tau2)
        if np.any(m2):
            S[m2] = np.exp(-((self.tau1 / self.eta0) ** self.beta0)) * np.exp(
                -self.lam1 * (t[m2] - self.tau1)
            )

        # Region 3: t >= tau2
        m3 = t >= self.tau2
        if np.any(m3):
            S[m3] = self.S_tau2 * np.exp(-(((t[m3] - self.tau2) / self.eta2) ** self.beta2))

        # t<0 => S=1 by convention
        return S

    def cond_fail_prob(self, t, dt):
        """P(fail in [t, t+dt) | survived to t). Works with scalars or arrays."""
        t = np.asarray(t, dtype=float)
        S_t = self.survival(t)
        S_tdt = self.survival(t + dt)
        p = 1.0 - np.clip(S_tdt / np.maximum(S_t, 1e-15), 0.0, 1.0)
        return p


def compute_fail_probs(phase, x_phase, n_sat, config):
    fail_model = config["fail_model"]
    launch_years = config["launch_years"]
    eval_years = config["eval_years"]

    fail_probs = np.zeros(n_sat)
    eval_year = eval_years[phase]

    for phase_past in range(phase + 1):
        # identify satellites launched in phase past
        sats_phase_past = np.where(x_phase == phase_past)[0]
        fail_probs[sats_phase_past] = 1 - fail_model.survival(eval_year - launch_years[phase_past])

    return fail_probs
