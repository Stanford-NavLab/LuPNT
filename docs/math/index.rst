Math Specifications
===================================================================

This section records the mathematical contracts used by LuPNT models.  The
goal is to make each model's state conventions, units, time scales, optional
corrections, and observable equations explicit enough that implementations can
be tested against the same specification.  Each page ties its algorithms to
the implementing C++/Python code.

These same contracts are honored by pure-Python model subclasses: Applications,
Agents, and Measurements can be authored in Python (pybind11 trampolines) and
registered into a YAML-driven ``Simulation`` via ``pnt.register_application`` /
``pnt.register_agent``, as demonstrated by Example 17.

.. toctree::
   :maxdepth: 1
   :caption: Fundamentals

   time_conversions
   frame_conversions
   dynamics
   integration
   autodiff

.. toctree::
   :maxdepth: 1
   :caption: Estimation

   filters

.. toctree::
   :maxdepth: 1
   :caption: Measurements & Signals

   gnss_measurements
   link_budget
   ionosphere_plasmasphere
   measurements

.. toctree::
   :maxdepth: 1
   :caption: Navigation Message Design

   ephemeris_almanac
