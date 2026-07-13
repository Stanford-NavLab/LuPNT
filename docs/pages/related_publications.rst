.. _related_publications:

Related Publications
====================

LuPNT is the simulation backbone for the Stanford NAV Lab's lunar Positioning,
Navigation, and Timing (PNT) research.  Some of the tutorial examples are related to the studies
published in the papers below.  This page maps each example (and the underlying
models) to the peer-reviewed work it is based on, so you can go from a notebook
to the methodology and results behind it.

Cite the simulator itself as:

.. code-block:: bibtex

   @inproceedings{IiyamaCasadesus2023,
     title     = {LuPNT: Open-Source Simulator for Lunar Positioning, Navigation, and Timing},
     author    = {Iiyama, Keidai and Casadesus Vila, Guillem and Gao, Grace},
     booktitle = {Proceedings of the Institute of Navigation GNSS+ conference (ION GNSS+ 2023)},
     year      = {2023},
     url       = {https://github.com/Stanford-NavLab/LuPNT},
   }

Examples mapped to publications
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Example(s)
     - Related publication(s)
   * - :doc:`ex1 <../tutorial/Python/ex1_propagate_orbit>`,
       :doc:`ex2 <../tutorial/Python/ex2_time_conversions>`
       (fundamentals: orbits, frames, time)
     - The dynamics, frame, and (lunicentric) time-scale conventions follow the
       author's PhD thesis [21]_ (Ch. 2) and the simulator papers [1]_ [13]_.
   * - :doc:`ex3 <../tutorial/Python/ex3_gnss_interface>`,
       :doc:`ex5 <../tutorial/Python/ex5_gnss_measurement_sim>`,
       :doc:`ex6 <../tutorial/Python/ex6_gnss_odts>`
       (terrestrial-GNSS interface, measurements, ODTS)
     - Lunar orbit & clock estimation from terrestrial GPS: [5]_ [18]_ (see
       :doc:`../math/filters`, :doc:`../math/gnss_measurements`).
   * - :doc:`ex4 <../tutorial/Python/ex4_plasmasphere>`
       (plasmasphere / ionosphere delay & ray tracing)
     - Iiyama & Gao, ionospheric and plasmaspheric delay characterization with
       GCPM: [20]_ [16]_ (arXiv:2510.10059).  See
       :doc:`../math/ionosphere_plasmasphere`.
   * - :doc:`ex7 <../tutorial/Python/ex7_groundstation_odts>`,
       :doc:`ex8 <../tutorial/Python/ex8_isl_odts>`
       (ground-station / inter-satellite-link ODTS)
     - Distributed lunar ODTS and time synchronization: [2]_ [7]_; ODTS with
       lunar surface stations: [15]_; the ISL-based autonomous clock-fault
       monitoring in ex8 is [19]_.
   * - :doc:`ex9 <../tutorial/Python/ex9_ephemeris>`
       (ephemeris & almanac fitting)
     - Ephemeris/almanac design for lunar navigation satellites: [9]_; earlier
       parameterization studies [6]_ [3]_ (see
       :doc:`../math/ephemeris_almanac`).
   * - :doc:`ex10 <../tutorial/Python/ex10_surface_rover>`
       (surface-rover / lander navigation)
     - Lunar surface user positioning: [4]_ [11]_; the full surface autonomy
       stack is [12]_ [17]_.
   * - :doc:`ex12 <../tutorial/Python/ex12_constellation_design>`
       (constellation design)
     - LANS constellation trade-off and staged deployment: [10]_.
   * - :doc:`ex15 <../tutorial/Python/ex15_marspnt>`
       (Mars PNT)
     - Orbit determination and time synchronization for a future Mars relay &
       navigation constellation: [14]_.
   * - General motivation / overview
     - "Can satellite-based radionavigation be extended to the Moon and other
       extraterrestrial bodies?" [8]_.

.. note::

   ``ex13`` (Cesium visualization) and ``ex16`` (LEO PNT with Harris-Priester
   drag) are simulator/demonstration examples without a dedicated publication.

Selected NAV Lab lunar-PNT publications
---------------------------------------

Ordered from oldest to newest.  A complete, up-to-date list is maintained on the
`Stanford NAV Lab publications page
<https://navlab.stanford.edu/publications/conference-articles>`_.

**2023**

.. [1] K. Iiyama, G. Casadesus Vila, and G. Gao, "`LuPNT: Open-source
   simulator for lunar positioning, navigation, and timing
   <https://doi.org/10.33012/2023.19373>`_," *ION GNSS+ 2023*.

.. [2] K. Iiyama and G. Gao, "`Positioning and timing of distributed lunar
   satellites via terrestrial GPS differential carrier phase measurements
   <https://doi.org/10.33012/2023.19385>`_," *ION GNSS+ 2023*.

.. [3] M. Cortinovis, K. Iiyama, and G. Gao, "`Satellite ephemeris
   approximation methods to support lunar positioning, navigation, and timing
   services <https://doi.org/10.33012/2023.19282>`_," *ION GNSS+ 2023* (Best
   Presentation of the Session).

.. [4] K. Iiyama, S. Bhamidipati, and G. Gao, "`Terrestrial GPS
   time-differenced carrier-phase positioning of lunar surface users
   <https://doi.org/10.1109/AERO55745.2023.10115673>`_," *2023 IEEE Aerospace
   Conference*.

**2024**

.. [5] K. Iiyama, S. Bhamidipati, and G. Gao, "`Precise positioning and
   timekeeping in a lunar orbit via terrestrial GPS time-differenced
   carrier-phase measurements <https://doi.org/10.33012/navi.635>`_,"
   *NAVIGATION: Journal of the Institute of Navigation*, 71(1), 2024.

.. [6] M. Cortinovis, K. Iiyama, and G. Gao, "`Satellite ephemeris
   parameterization methods to support lunar positioning, navigation, and
   timing services <https://doi.org/10.33012/navi.664>`_," *NAVIGATION: Journal
   of the Institute of Navigation*, 71(4), 2024.

.. [7] K. Iiyama, G. Casadesus Vila, and G. Gao, "`Contact plan optimization
   and distributed state estimation for delay tolerant satellite networks
   <https://doi.org/10.1109/AERO58975.2024.10521114>`_," *2024 IEEE Aerospace
   Conference*.

**2025**

.. [8] K. Iiyama, S. Pullen, and G. Gao, "`Can satellite-based radionavigation
   be extended to the Moon and other extraterrestrial bodies?
   <https://insidegnss.com/can-satellite-based-radionavigation-be-extended-to-the-moon-and-other-extraterrestrial-bodies/>`_"
   *Inside GNSS*, 20(5), 18-27, 2025.

.. [9] K. Iiyama and G. Gao, "`Ephemeris and almanac design for lunar
   navigation satellites <https://arxiv.org/abs/2510.25161>`_," *IEEE
   Transactions on Aerospace and Electronic Systems*, under review, 2025.

.. [10] K. Iiyama and G. Gao, "`Trade-off analysis for lunar augmented
   navigation service constellation design <https://arxiv.org/abs/2510.16030>`_,"
   *NAVIGATION: Journal of the Institute of Navigation*, under review, 2025.

.. [11] K. M. Y. Coimbra, M. Cortinovis, T. Mina, and G. Gao, "`Single-satellite
   lunar navigation via Doppler shift observables for the NASA Endurance mission
   <https://doi.org/10.33012/navi.710>`_," *NAVIGATION: Journal of the Institute
   of Navigation*, 72(3), 2025.

.. [12] A. Dai, G. Casadesus Vila, A. Wu, K. Iiyama, K. Coimbra, T. Deng, and
   G. Gao, "Full stack navigation, mapping, and planning for the lunar autonomy
   challenge," *NAVIGATION: Journal of the Institute of Navigation*, under
   review, 2025 (conference version: [17]_).

.. [13] G. Casadesus Vila, K. Iiyama, and G. Gao, "`LuPNT: An open-source
   simulator for lunar communications, positioning, navigation, and timing
   <https://doi.org/10.1109/AERO63441.2025.11068501>`_," *2025 IEEE Aerospace
   Conference*.

.. [14] K. Iiyama, W. W. Jun, S. Bhamidipati, G. Gao, and K. Cheung, "`Orbit
   determination and time synchronization for the future Mars relay and
   navigation constellation <https://doi.org/10.1109/AERO63441.2025.11068793>`_,"
   *2025 IEEE Aerospace Conference*.

.. [15] G. Casadesus Vila and G. Gao, "`Moon surface station to support lunar
   positioning, navigation, and timing services
   <https://doi.org/10.33012/2025.20258>`_," *ION GNSS+ 2025*.

.. [16] K. Iiyama and G. Gao, "`Ionospheric and plasmaspheric delay
   characterization and mitigation methodologies for lunar terrestrial GNSS
   receivers <https://doi.org/10.33012/2025.20343>`_," *ION GNSS+ 2025* (Best
   Presentation of the Session).

.. [17] A. Dai, A. Wu, K. Iiyama, G. Casadesus Vila, K. Coimbra, T. Deng, and
   G. Gao, "`Full stack navigation, mapping, and planning for the lunar autonomy
   challenge <https://doi.org/10.33012/2025.20447>`_," *ION GNSS+ 2025*.

**2026**

.. [18] K. Iiyama and G. Gao, "`GNSS-based lunar orbit and clock estimation with
   stochastic cloning UD filter <https://arxiv.org/abs/2601.16393>`_," *Journal
   of Guidance, Control, and Dynamics*, under review, 2026.

.. [19] K. Iiyama, D. Neamati, and G. Gao, "`Satellite autonomous clock fault
   monitoring with inter-satellite ranges using Euclidean distance matrices
   <https://doi.org/10.33012/navi.764>`_," *NAVIGATION: Journal of the Institute
   of Navigation*, 73(1), 2026.

.. [20] K. Iiyama and G. Gao, "`Ionospheric and plasmaspheric delay
   characterization for lunar terrestrial GNSS receivers with Global Core Plasma
   Model <https://arxiv.org/abs/2510.10059>`_," *NAVIGATION: Journal of the
   Institute of Navigation*, accepted, 2026.

.. [21] K. Iiyama, "`Design and algorithms for lunar navigation satellite
   systems <https://purl.stanford.edu/wx922dp2089>`_," *Ph.D. Thesis, Stanford
   University*, 2026.
