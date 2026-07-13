.. _cpp_gnss_tutorials:

GNSS: Interface, Measurements, and ODTS
=======================================

Three of the example programs form the Earth-GNSS chain — the constellation
interface, receiver-side measurement simulation, and lunar orbit determination.
Each mirrors the correspondingly numbered Python notebook and has its own
tutorial page (with the full C++ source embedded):

* :doc:`Example 3 — GNSS Interface <ex3_gnss_interface>` — constellation catalog,
  transmit-antenna gain patterns, and a broadcast-vs-precise ephemeris comparison
  (SP3 vs. BRDC with phase-centre-offset corrections).
* :doc:`Example 5 — GNSS Measurement Simulation <ex5_gnss_measurement_sim>` —
  propagate an ELFO receiver and run visibility + link-budget (C/N\ :sub:`0`)
  against the GPS + Galileo constellations.
* :doc:`Example 6 — GNSS Sidelobe ODTS <ex6_gnss_odts>` — sidelobe pseudorange +
  TDCP orbit determination and time synchronization with a UDU EKF, including a
  GCPM/IRI plasmaspheric-delay ray-trace stage.

Shared data setup
-----------------

All three are driven by SP3 precise-ephemeris products (and ANTEX / BRDC), which
require a free NASA Earthdata login and a one-time download — see
:doc:`../../pages/sp3_download`. ``pixi run build-examples`` only *compiles* the
programs; the data is loaded at run time. Each program is written to guard
missing data gracefully — it prints an informative note and exits cleanly when
the SP3 / BRDC / ANTEX inputs (or Earthdata credentials) are absent — so it
always builds and runs, and produces its full analysis once the data is present.
Example 6 additionally needs the staged link + plasma-delay precompute produced
by ``python/examples/ex6_precompute.py``.

The detailed receiver-side measurement model shared by ex5 and ex6 (light-time
iteration, relativistic and plasma corrections, antenna gain, and the online vs.
batch ``Precompute`` paths) is documented on the :doc:`gnss_measurements` page.
