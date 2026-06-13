import numpy as np


class ClockNoise:
    def __init__(self, clock_model):
        self.clock_model = clock_model
        if clock_model == "OCXO":
            # For LCRNS
            self.sigma_wf = 6.2299445014e-13  # White frequency noise (s/sqrt(s))
            self.sigma_rw = 2.0129544799e-14  # Random walk frequency noise (s^(-1/2))
            self.sigma_rr = 7.0118586804e-28  # Random run frequency noise (s^(-3/2))
            self.sigma_wp = 1.4562174977e-19  # White phase noise (s)
            self.sigma_z0 = 3.1943096158e-16  # Initial clock bias (s^(-1))
        else:
            raise ValueError(f"Unsupported clock model: {clock_model}")

    def get_process_noise(self, tau):
        s2_wf = self.sigma_wf**2
        s2_rw = self.sigma_rw**2
        s2_rr = self.sigma_rr**2

        q11 = s2_wf * tau + (s2_rw * tau**3) / 3 + (s2_rr * tau**5) / 20
        q12 = (s2_rw * tau**2) / 2 + (s2_rr * tau**4) / 8
        q13 = (s2_rr * tau**3) / 6
        q22 = s2_rw * tau + (s2_rr * tau**3) / 3
        q23 = (s2_rr * tau**2) / 2
        q33 = s2_rr * tau
        Q = np.array([[q11, q12, q13], [q12, q22, q23], [q13, q23, q33]])
        return Q

    def simulate_clock_bias(self, tspan):
        n = len(tspan)
        dt = np.diff(tspan, prepend=0)
        clk_bias = np.zeros((n, 3))  # [bias, drift, drift_rate]

        # sample initial clock bias
        clk_bias[0, 0] = np.random.normal(0, self.sigma_z0)

        for i in range(1, n):
            tau = dt[i]
            Q = self.get_process_noise(tau)
            w = np.random.multivariate_normal(mean=[0, 0, 0], cov=Q)
            Phi = np.array([[1, tau, 0.5 * tau**2], [0, 1, tau], [0, 0, 1]])
            clk_bias[i] = Phi @ clk_bias[i - 1] + w

        return clk_bias
