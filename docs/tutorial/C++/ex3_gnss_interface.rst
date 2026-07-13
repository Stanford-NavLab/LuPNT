.. _cpp_ex3_gnss_interface:

Example 3: GNSS Interface
=========================

LuPNT's Earth-GNSS data interfaces: load the GPS / Galileo / QZSS constellations,
inspect transmit-antenna gain patterns (ANTEX / ACE), and compare broadcast
versus precise ephemerides (SP3 vs. BRDC, with phase-centre-offset corrections
and a radial / along-track / cross-track decomposition of the difference).

Mirrors :doc:`../Python/ex3_gnss_interface`.

In C++ the constellation catalog is ``GnssConstellation`` (``lupnt/agents``),
antenna patterns are loaded by ``Antenna`` and phase-centre offsets by
``AntexLoader`` (``lupnt/interfaces``), and the precise / broadcast ephemerides
are parsed by ``Sp3Loader`` / ``RinexNavLoader``. The program prints the
transmit-antenna peak gains, the ANTEX phase-centre offsets, and a per-PRN
radial / along-track / cross-track + clock error summary of BRDC against
(SP3 + PCO).

.. note::

   This example is data-driven: it reads SP3 precise-ephemeris products (plus
   ANTEX and BRDC), which require a free NASA Earthdata login and a one-time
   download (see :doc:`../../pages/sp3_download`). ``pixi run build-examples``
   only *compiles* the program; the data is loaded at run time, and the program
   guards missing inputs — it prints an informative note and exits cleanly when
   the SP3 / BRDC / ANTEX files or Earthdata credentials are absent, so it always
   builds and runs and produces its full analysis once the data is present.

.. literalinclude:: ../../../cpp/examples/tutorials/ex3_gnss_interface.cc
   :language: cpp

.. code-block:: bash

   pixi run build-examples
   ./build-examples-pixi/ex3_gnss_interface
