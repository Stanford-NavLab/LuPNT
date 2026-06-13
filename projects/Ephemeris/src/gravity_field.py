import numpy as np
import math
import re


class GravityField:
    def __init__(self, n, m):
        self.n_max = 0
        self.m_max = 0
        self.n = n
        self.m = m
        self.GM = 0
        self.R = 0
        self.CS = np.zeros((n + 1, m + 1))
        self.CS[0, 0] = 1.0


def kron(i, j):
    return 1 if i == j else 0


def factprod(n, m):
    result = math.factorial(n - m) / math.factorial(n + m)
    return result


KM_M = 1e-3  # Conversion from km to meters


def read_harmonic_gravity_field(filename, n, m, normalized):
    gravity_field = GravityField(n, m)

    with open(filename, "r") as file:
        for line in file:
            if line.startswith("POTFIELD"):
                n_max_in = int(line[8:11].strip())
                m_max_in = int(line[11:14].strip())
                parts = line.split()
                GM = float(parts[-3])
                r = float(parts[-2])

                gravity_field.n_max = n_max_in
                gravity_field.m_max = m_max_in
                gravity_field.GM = GM * KM_M**3  # KM to meters
                gravity_field.R = r * KM_M
                break

        for line in file:
            if line.startswith("RECOEF"):
                line_content = line[6:].strip()
                n_in = int(line_content[:3].strip())
                m_in = int(line_content[3:6].strip())

                if n_in > n:
                    continue
                if m_in > m:
                    continue

                coef_content = line_content[6:].strip()
                coef_parts = re.findall(r"[+-]?\d+\.\d+[eE][+-]?\d+", coef_content)

                cnm = float(coef_parts[0])
                snm = float(coef_parts[1]) if len(coef_parts) > 1 else 0.0

                if m_in == 0:
                    N = math.sqrt(2 * n_in + 1) if normalized else 1.0
                    gravity_field.CS[n_in, m_in] = N * cnm
                else:
                    N = (
                        math.sqrt((2 - kron(0, m_in)) * (2 * n_in + 1) * factprod(n_in, m_in))
                        if normalized
                        else 1.0
                    )
                    gravity_field.CS[n_in, m_in] = N * cnm
                    gravity_field.CS[m_in - 1, n_in] = N * snm

    return gravity_field


def get_cnm(CS, n, m):
    return CS[n, m]


def get_snm(CS, n, m):
    if m == 0:
        # error
        return 0
    return CS[m - 1, n]


def get_cnm_mat(CS, n, m):
    CS_mat = np.zeros(CS.shape)
    CS_mat[n, m] = CS[n, m]
    return CS_mat


def get_snm_mat(CS, n, m):
    if m == 0:
        # error
        return 0
    CS_mat = np.zeros(CS.shape)
    CS_mat[m - 1, n] = CS[m - 1, n]
    return CS_mat
