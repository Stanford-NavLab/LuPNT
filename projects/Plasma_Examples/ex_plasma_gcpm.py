"""
ex_plasma_gcpm.py — GCPM v2.4 example (Python script version of ex_plasma_gcpm.ipynb)

Run inside the pixi environment:
    pixi run python projects/Plasma_Examples/ex_plasma_gcpm.py
"""

import numpy as np
from pylupnt import plasma as pecsimpy
from pylupnt.plasma.kp_loader import update_kp

# ── Setup ────────────────────────────────────────────────────────────────────
base_path = pecsimpy.get_plasma_base_path()
print("Base path:", base_path)
update_kp(base_path)

year, doy, hour, minute, sec = 2002, 185, 12, 0, 0
dt_cpp = pecsimpy.DateTime(year, doy, hour, minute, sec)

r = 1.014
akp = 0.7
al = 1.107226364
alatr = np.arccos(np.sqrt(r / al))  # magnetic latitude [rad]
amlt = 23.74  # magnetic local time [hours]

pecsimpy.set_iri_model("IRI2007")

mjd = pecsimpy.datetime_to_mjd(dt_cpp)
tj2000 = pecsimpy.mjd_to_tj2000(mjd)

iri2007option = pecsimpy.IRI2007Option()
iri2007option.R12 = -1.0  # -1 → use IRI internally
pecsimpy.set_iri2007_option(iri2007option)

ionoparams = pecsimpy.get_iono_params(tj2000, akp)

# ── Case 1 ───────────────────────────────────────────────────────────────────
out_cpp = pecsimpy.gcpm_v24(dt_cpp, r, amlt, alatr, akp)
out_fortran = pecsimpy.gcpm_v24_fortran(dt_cpp, r, amlt, alatr, akp)
print("\n[Case 1] Default case")
print(f"  Year={year}, DOY={doy}, Hour={hour}, Min={minute}, Sec={sec}")
print(f"  alatr={alatr:.4f} rad, amlt={amlt}, r={r} RE")
print(f"  F10.7={out_cpp[4]:.2f}  (ionoparams={ionoparams[0]:.2f})")
print(f"  R12={out_cpp[5]:.2f}    (ionoparams={ionoparams[1]:.2f})")
print(f"  Ne, H+, He+, O+ [cm^-3] (C++):     {out_cpp[:4]}")
print(f"  Ne, H+, He+, O+ [cm^-3] (Fortran):  {out_fortran}")

# ── Case 2 ───────────────────────────────────────────────────────────────────
r += 0.5
out_cpp = pecsimpy.gcpm_v24(dt_cpp, r, amlt, alatr, akp)
out_fortran = pecsimpy.gcpm_v24_fortran(dt_cpp, r, amlt, alatr, akp)
print("\n[Case 2] Larger radius")
print(f"  r={r} RE")
print(f"  Ne, H+, He+, O+ [cm^-3] (C++):     {out_cpp[:4]}")
print(f"  Ne, H+, He+, O+ [cm^-3] (Fortran):  {out_fortran}")

# ── Case 3 ───────────────────────────────────────────────────────────────────
dt_cpp.year, dt_cpp.doy, dt_cpp.hour, dt_cpp.min, dt_cpp.sec = 2025, 1, 12, 0, 0
r, akp = 2.2, 6.5
out_cpp = pecsimpy.gcpm_v24(dt_cpp, r, amlt, alatr, akp)
out_fortran = pecsimpy.gcpm_v24_fortran(dt_cpp, r, amlt, alatr, akp)
print("\n[Case 3] Different date / radius / Kp")
print(f"  Year={dt_cpp.year}, DOY={dt_cpp.doy}, r={r} RE, akp={akp}")
print(f"  Ne, H+, He+, O+ [cm^-3] (C++):     {out_cpp[:4]}")
print(f"  Ne, H+, He+, O+ [cm^-3] (Fortran):  {out_fortran}")

# ── Case 4 ───────────────────────────────────────────────────────────────────
akp = pecsimpy.get_kp_index(dt_cpp)
out_cpp = pecsimpy.gcpm_v24(dt_cpp, r, amlt, alatr, akp)
out_fortran = pecsimpy.gcpm_v24_fortran(dt_cpp, r, amlt, alatr, akp)
print("\n[Case 4] Kp from index file")
print(f"  akp={akp:.2f}")
print(f"  Ne, H+, He+, O+ [cm^-3] (C++):     {out_cpp[:4]}")
print(f"  Ne, H+, He+, O+ [cm^-3] (Fortran):  {out_fortran}")

# ── Case 5 ───────────────────────────────────────────────────────────────────
iri2007option.R12 = 1.0
pecsimpy.set_iri2007_option(iri2007option)
out_cpp = pecsimpy.gcpm_v24(dt_cpp, r, amlt, alatr, akp)
out_fortran = pecsimpy.gcpm_v24_fortran(dt_cpp, r, amlt, alatr, akp)
print("\n[Case 5] R12=1.0 override")
print(f"  F10.7={out_cpp[4]:.2f}, R12={out_cpp[5]:.2f}")
print(f"  Ne, H+, He+, O+ [cm^-3] (C++):     {out_cpp[:4]}")
print(f"  Ne, H+, He+, O+ [cm^-3] (Fortran):  {out_fortran}")

# ── Case 6 ───────────────────────────────────────────────────────────────────
pecsimpy.set_iri_model("IRI2020")
out_cpp = pecsimpy.gcpm_v24(dt_cpp, r, amlt, alatr, akp)
out_fortran = pecsimpy.gcpm_v24_fortran(dt_cpp, r, amlt, alatr, akp)
print("\n[Case 6] IRI2020 model")
print(f"  Ne, H+, He+, O+ [cm^-3] (C++ + IRI2020):    {out_cpp[:4]}")
print(f"  Ne, H+, He+, O+ [cm^-3] (Fortran + IRI2007): {out_fortran}")
