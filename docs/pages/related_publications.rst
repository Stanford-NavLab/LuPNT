.. _related_publications:

Related Publications
====================

LuPNT is the simulation backbone for the Stanford NAV Lab's lunar Positioning,
Navigation, and Timing (PNT) research, most of it by Keidai Iiyama and
collaborators.  The tutorial examples reproduce, in miniature, the studies
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
       author's PhD thesis, Ch. 2, and the simulator papers [C16]_ [C20]_.
   * - :doc:`ex3 <../tutorial/Python/ex3_gnss_interface>`,
       :doc:`ex5 <../tutorial/Python/ex5_gnss_measurement_sim>`,
       :doc:`ex6 <../tutorial/Python/ex6_gnss_odts>`
       (terrestrial-GNSS interface, measurements, ODTS)
     - Lunar orbit & clock estimation from terrestrial GPS: [J2]_ [C15]_ [C11]_;
       the stochastic-cloning UD filter used for TDCP ODTS is [J9]_ (see
       :doc:`../math/filters`, :doc:`../math/gnss_measurements`).
   * - :doc:`ex4 <../tutorial/Python/ex4_plasmasphere>`
       (plasmasphere / ionosphere delay & ray tracing)
     - Iiyama & Gao, ionospheric and plasmaspheric delay characterization with
       GCPM: [J6]_ [C22]_ (arXiv:2510.10059).  See
       :doc:`../math/ionosphere_plasmasphere`.
   * - :doc:`ex7 <../tutorial/Python/ex7_groundstation_odts>`,
       :doc:`ex8 <../tutorial/Python/ex8_isl_odts>`
       (ground-station / inter-satellite-link ODTS)
     - Distributed lunar ODTS and time synchronization: [C7]_; the ISL-based
       autonomous clock-fault monitoring in ex8 is [J4]_ [C18]_.
   * - :doc:`ex9 <../tutorial/Python/ex9_ephemeris>`
       (ephemeris & almanac fitting)
     - Ephemeris/almanac design for lunar navigation satellites: [J7]_; earlier
       parameterization studies [J3]_ [C14]_ (see
       :doc:`../math/ephemeris_almanac`).
   * - :doc:`ex10 <../tutorial/Python/ex10_surface_rover>`,
       :doc:`ex11 <../tutorial/Python/ex11_lander_navigation>`
       (surface-rover / lander navigation)
     - Lunar surface user positioning and time transfer: [C12]_ [C9]_; the full
       surface autonomy stack is [J8]_ [C23]_.
   * - :doc:`ex12 <../tutorial/Python/ex12_constellation_design>`
       (constellation design)
     - LANS constellation trade-off and staged deployment: [J5]_ [C21]_; early
       constellation optimization [C2]_.
   * - :doc:`ex14 <../tutorial/Python/ex14_opnav>`
       (optical navigation)
     - Autonomous angles-only navigation and timekeeping in lunar orbit: [C8]_.
   * - :doc:`ex15 <../tutorial/Python/ex15_marspnt>`
       (Mars PNT)
     - Orbit determination and time synchronization for a future Mars relay &
       navigation constellation: [C19]_.
   * - General motivation / overview
     - "Can satellite-based radionavigation be extended to the Moon and other
       extraterrestrial bodies?" [M1]_.

.. note::

   ``ex13`` (Cesium visualization) and ``ex16`` (LEO PNT with Harris-Priester
   drag) are simulator/demonstration examples without a dedicated publication.

Selected NAV Lab lunar-PNT publications
---------------------------------------

.. [M1] K. Iiyama, S. Pullen, and G. Gao, "Can satellite-based radionavigation
   be extended to the Moon and other extraterrestrial bodies?" *Inside GNSS*,
   20(5), 18-27, 2025.

.. [J9] K. Iiyama and G. Gao, "GNSS-based lunar orbit and clock estimation with
   stochastic cloning UD filter," *Journal of Guidance, Control, and Dynamics*,
   under review, 2026.

.. [J7] K. Iiyama and G. Gao, "Ephemeris and almanac design for lunar
   navigation satellites," *IEEE Transactions on Aerospace and Electronic
   Systems*, under review, 2025.

.. [J6] K. Iiyama and G. Gao, "Ionospheric and plasmaspheric delay
   characterization for lunar terrestrial GNSS receivers with Global Core
   Plasma Model," *NAVIGATION: Journal of the Institute of Navigation*, 2025.
   arXiv:2510.10059.

.. [J5] K. Iiyama and G. Gao, "Trade-off analysis for lunar augmented
   navigation service constellation design," *NAVIGATION: Journal of the
   Institute of Navigation*, under review, 2025.

.. [J4] K. Iiyama, D. Neamati, and G. Gao, "Satellite autonomous clock fault
   monitoring with inter-satellite ranges using Euclidean distance matrices,"
   *NAVIGATION: Journal of the Institute of Navigation*, 73(1), 2026.

.. [J3] M. Cortinovis, K. Iiyama, and G. Gao, "Satellite ephemeris
   parameterization methods to support lunar positioning, navigation, and
   timing services," *NAVIGATION: Journal of the Institute of Navigation*,
   71(4), 2024.

.. [J2] K. Iiyama, S. Bhamidipati, and G. Gao, "Precise positioning and
   timekeeping in a lunar orbit via terrestrial GPS time-differenced
   carrier-phase measurements," *NAVIGATION: Journal of the Institute of
   Navigation*, 71(1), 2024.

.. [J8] A. Dai, G. Casadesus Vila, A. Wu, K. Iiyama, K. Coimbra, T. Deng, and
   G. Gao, "Full stack navigation, mapping, and planning for the lunar autonomy
   challenge," *NAVIGATION: Journal of the Institute of Navigation*, under
   review, 2025.

.. [C23] A. Dai, A. Wu, K. Iiyama, G. Casadesus Vila, K. Coimbra, T. Deng, and
   G. Gao, "Full stack navigation, mapping, and planning for the lunar autonomy
   challenge," *ION GNSS+ 2025*.

.. [C22] K. Iiyama and G. Gao, "Ionospheric and plasmaspheric delay
   characterization and mitigation methodologies for lunar terrestrial GNSS
   receivers," *ION GNSS+ 2025* (Best Presentation of the Session).

.. [C21] K. Iiyama and G. Gao, "Constellation design and staged development for
   the lunar navigation satellite system," *ION GNSS+ 2025*.

.. [C20] G. Casadesus Vila, K. Iiyama, and G. Gao, "LuPNT: An open-source
   simulator for lunar communications, positioning, navigation, and timing,"
   *2025 IEEE Aerospace Conference*.

.. [C19] K. Iiyama, W. W. Jun, S. Bhamidipati, G. Gao, and K. Cheung, "Orbit
   determination and time synchronization for the future Mars relay and
   navigation constellation," *2025 IEEE Aerospace Conference*.

.. [C18] K. Iiyama, D. Neamati, and G. Gao, "Autonomous constellation fault
   monitoring with inter-satellite links: A rigidity-based approach,"
   *ION GNSS+ 2024* (Best Presentation of the Session).

.. [C16] K. Iiyama, G. Casadesus Vila, and G. Gao, "LuPNT: Open-source
   simulator for lunar positioning, navigation, and timing," *ION GNSS+ 2023*.

.. [C15] K. Iiyama and G. Gao, "Positioning and timing of distributed lunar
   satellites via terrestrial GPS differential carrier phase measurements,"
   *ION GNSS+ 2023*.

.. [C14] M. Cortinovis, K. Iiyama, and G. Gao, "Satellite ephemeris
   approximation methods to support lunar positioning, navigation, and timing
   services," *ION GNSS+ 2023* (Best Presentation of the Session).

.. [C12] K. Iiyama, S. Bhamidipati, and G. Gao, "Terrestrial GPS
   time-differenced carrier-phase positioning of lunar surface users,"
   *2023 IEEE Aerospace Conference*.

.. [C11] K. Iiyama, S. Bhamidipati, and G. Gao, "Precise positioning and
   timekeeping in a lunar orbit via terrestrial GPS time-differenced
   carrier-phase measurements," *2023 International Technical Meeting of the
   Institute of Navigation*.

.. [C9] S. Bhamidipati, K. Iiyama, T. Mina, and G. Gao, "Time-transfer from
   terrestrial GPS for distributed lunar surface communication networks,"
   *2022 IEEE Aerospace Conference*.

.. [C8] K. Iiyama, J. Kruger, and S. D'Amico, "Autonomous distributed
   angles-only navigation and timekeeping in lunar orbit," *2022 International
   Technical Meeting of the Institute of Navigation*.

.. [C7] K. Iiyama, Y. Kawabata, and R. Funase, "Autonomous and decentralized
   orbit determination and clock offset estimation of lunar navigation
   satellites using GPS signals and inter-satellite ranging," *ION GNSS+ 2021*.

.. [C2] K. Iiyama, "Optimization of the navigation satellite constellation and
   lunar monitoring station for lunar global navigation satellite system,"
   *32nd International Symposium on Space Technology and Science*, 2019.

A complete, up-to-date list is maintained at
`kdricemt.github.io/publications <https://kdricemt.github.io/publications/>`_
and on the `Stanford NAV Lab publications page
<https://navlab.stanford.edu/publications/conference-articles>`_.
