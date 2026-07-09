Automatic Differentiation
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines how LuPNT computes derivatives -- state-transition
matrices, measurement Jacobians, and force-model sensitivities -- by
**automatic differentiation (AD)** rather than by hand-coded partials or finite
differences.  AD is the mechanism behind every :math:`\Phi = \partial x_f/
\partial x_0` in :doc:`dynamics` / :doc:`integration`, every
:math:`H = \partial z/\partial x` in :doc:`gnss_measurements` /
:doc:`measurements`, and the linearized dynamics/measurement callbacks
consumed by the :doc:`filters`.  The goal is to make the differentiable-type
contract, the Jacobian API, and the seeding conventions explicit enough that
new models get exact derivatives "for free" and tests can compare against the
analytic expressions.

The library and the differentiable scalar are declared in
``cpp/lupnt/core/definitions.h``.

The Differentiable Scalar ``Real``
-------------------------------------------------------------------

LuPNT builds on the `autodiff <https://autodiff.github.io>`_ library in
**forward mode**.  The fundamental scalar is

.. math::

   \texttt{Real} \equiv \texttt{autodiff::real},

a first-order forward-mode dual number carrying a value and one directional
derivative, :math:`(\,\texttt{val},\ \texttt{grad}\,)`.  Every LuPNT vector and
matrix type is an Eigen container over ``Real`` (``Vec3``, ``Vec6``, ``VecX``,
``Mat3``, ``MatX``, ...), with plain-``double`` counterparts (``Vec3d``,
``MatXd``, ...) for storage and I/O:

.. code-block:: cpp

   // cpp/lupnt/core/definitions.h
   namespace ad = autodiff;
   using Real = ad::real;                 // forward-mode, first order
   using ad::at;  using ad::wrt;  using ad::jacobian;  using ad::derivative;
   template <int size> using Vec = Matrix<Real, size, 1>;   // Eigen over Real
   template <int r, int c> using Mat = Matrix<Real, r, c>;

Because the models (dynamics ``ComputeRates``, force kernels in
``environment/forces``, measurement ``Compute``, frame rotations, ...) are all
written in terms of ``Real`` arithmetic, a derivative propagates through any of
them without additional code.  The numeric value of a ``Real`` is recovered
with ``.val()``; ``cast<double>()`` converts a ``Real`` Eigen expression to its
``double`` values.

.. note::

   Forward mode is first-order only: it yields Jacobians, not Hessians, and the
   propagated state-transition matrix is therefore the first-order variational
   :math:`\Phi = \partial x_f/\partial x_0` (see :doc:`dynamics`).  Because
   ``autodiff::real`` carries a single derivative direction, a full :math:`m
   \times n` Jacobian costs :math:`n` forward sweeps (one per input), which
   ``jacobian`` and ``JacobianParallel`` (below) organize column by column.

The Jacobian Contract
-------------------------------------------------------------------

The single AD entry point is autodiff's ``jacobian``, which evaluates a
vector-valued function and its Jacobian in the same call.  For
:math:`y = f(x)` with :math:`x \in \mathbb{R}^n`, :math:`y \in \mathbb{R}^m`,

.. math::

   J = \frac{\partial y}{\partial x} \in \mathbb{R}^{m\times n},
   \qquad
   \texttt{jacobian}(f,\ \texttt{wrt}(x),\ \texttt{at}(x),\ y,\ J),

where ``wrt(x)`` names the differentiation variables (all of ``x``, or a single
component ``x(i)``, or a parameter block ``p``) and ``at(...)`` supplies the
full evaluation point (which may include additional non-differentiated
arguments such as :math:`t_0, t_f, u`).  The value :math:`y` and the Jacobian
:math:`J` are returned together.

A recurring idiom is to **reset the seeds** by copying the input through
``cast<double>()`` before seeding it, so that the Jacobian is taken cleanly with
respect to the requested variables and no stale derivative directions leak in:

.. code-block:: cpp

   VecX x_tmp = x.cast<double>();          // drop any incoming seeds
   jacobian(f, wrt(x_tmp), at(x_tmp), y, J);   // J = df/dx at x

State-Transition Matrix (Dynamics)
-------------------------------------------------------------------

The default ``Dynamics::Propagate(x0, t0, tf, u, stm)`` forms the STM by
differentiating the (seed-reset) no-STM propagation with respect to the initial
state (``cpp/lupnt/dynamics/dynamics.cc``):

.. code-block:: cpp

   VecX x0_tmp = x0.cast<double>();
   jacobian(func, wrt(x0_tmp), at(x0_tmp), xf_tmp, *stm);   // stm = d(xf)/d(x0)

where ``func`` re-propagates the state.  This is the :math:`\Phi(t_f,t_0)` of
:doc:`dynamics`; for ``NBodyDynamics`` the STM overload requires AD to be
enabled (``use_ad_``).  The joint state/parameter overload
``PropagateWithParams`` runs two sweeps to build the two sensitivity blocks
(used by batch/EKF estimation of force-model or clock parameters):

.. code-block:: cpp

   jacobian(func, wrt(x0_tmp), at(x0_tmp, p0_tmp), xf_tmp, *stm_state);  // d(xf)/d(x0)
   jacobian(func, wrt(p0_tmp), at(x0_tmp, p0_tmp), xf_tmp, *stm_param);  // d(xf)/d(params)

Integrator Sensitivity (column-parallel AD)
-------------------------------------------------------------------

For the integrator sensitivity matrix, ``JacobianParallel``
(``cpp/lupnt/numerics/math_utils.cc``) differentiates the re-propagation
function **one input at a time**, distributing the columns across cores with
OpenMP.  Each thread seeds a single component ``wrt(x0(i))`` and runs one
forward sweep to produce the corresponding column of :math:`J`:

.. code-block:: cpp

   #pragma omp parallel for
   for (int i = 0; i < n; i++) {
     VecX  xi;  MatXd Ji;
     VecX  x0_tmp = x0.cast<double>();
     jacobian(func, wrt(x0_tmp(i)), at(x0_tmp), xi, Ji);   // one column
     J_vec[i] = Ji;
   }

.. note::

   ``JacobianParallel`` is exact **forward-mode automatic differentiation**
   evaluated column-by-column, not a finite-difference approximation: ``func``
   propagates in ``Real`` arithmetic, so the derivative is carried through the
   integration analytically.  The cost is :math:`n` forward sweeps (re-runs of
   the trajectory), which is why the columns are parallelized.

Filter Dynamics and Measurement Jacobians
-------------------------------------------------------------------

Plain dynamics/measurement models (returning only the value) are lifted into
the filter callbacks -- which must also return the linearization :math:`F` and
:math:`H` -- by wrapping them in ``jacobian``
(``cpp/lupnt/numerics/filters/filter_utils.cc``):

.. code-block:: cpp

   // GetFilterMeasurementFunction: H = dz/dx
   FilterMeasurementFunction f = [f_meas](const State& x, MatXd* H, MatXd* R) -> VecX {
     VecX y;  VecX x_tmp = x.cast<double>();
     jacobian(f_meas, wrt(x_tmp), at(x_tmp, R), y, *H);
     return y;
   };

   // GetFilterDynamicsFunction: F = d(xf)/d(x0)
   FilterDynamicsFunction g = [f_dyn](const State& x, Real t0, Real tf, const State* u, MatXd* F) {
     VecX xf;  VecX x_tmp = x.cast<double>();
     jacobian(f_dyn, wrt(x_tmp), at(x_tmp, t0, tf, u), xf, *F);
     return State(xf, x.GetNames(), x.GetUnits());
   };

These are exactly the ``f_dyn_`` / ``f_meas_`` callbacks the :doc:`filters`
consume (:math:`F` in the EKF predict, :math:`H` in the update).  A measurement
model that instead differentiates its own raw ``VecX`` model in-place -- e.g.
``IslCrosslinkMeasurement`` (:doc:`measurements`) -- calls ``jacobian`` directly
and deliberately re-seeds a fresh state so that the incoming derivative seeds
are preserved for the ``SchmidtEKF`` consider-state update.

Differentiability Requirements for Models
-------------------------------------------------------------------

For AD to yield finite, correct derivatives, model code must stay
differentiable through every evaluated branch:

* Write model math in ``Real`` (not ``double``); use ``.val()`` only for
  logging, comparisons, or when returning to non-differentiated interfaces.
* Guard non-smooth boundaries.  At the domain edges of ``asin``/``acos`` the
  derivative is infinite, so the SRP shadow function clamps its arguments to a
  strictly interior value (``ClampUnit``, :math:`\pm(1-10^{-12})`) so that the
  SRP contribution to :math:`\Phi` stays finite across a grazing eclipse
  boundary (see the SRP note in :doc:`dynamics`).
* Avoid discontinuous control flow keyed on the differentiated variable where a
  smooth derivative is expected; ``jacobian`` returns the derivative of the
  branch actually taken.

Model Boundaries
-------------------------------------------------------------------

* AD is **forward mode, first order** (``autodiff::real``): Jacobians only, no
  second derivatives; a full :math:`m\times n` Jacobian costs :math:`n` forward
  sweeps.
* The ``cast<double>()`` seed reset is required before an independent
  ``jacobian`` call; omitting it can mix stale derivative directions into the
  result.
* ``NBodyDynamics::Propagate(..., stm)`` requires ``use_ad_`` enabled; models
  that are not written in ``Real`` cannot be differentiated this way.
* AD differentiates the code as written -- non-smooth clamps, table lookups, or
  branch selections produce the derivative of the realized path, which is why
  smoothing helpers such as ``ClampUnit`` are used at known singular edges.
