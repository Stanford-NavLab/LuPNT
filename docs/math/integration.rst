Numerical Integration
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT's ODE
integrators, which advance a state under ``dx/dt = f(t, x)`` and, when
requested, propagate the associated sensitivity (state-transition) matrix.
The goal is to make the step equations, embedded error estimates, adaptive
step-size control, and termination semantics explicit enough that propagation
and filtering tests can compare against the same equations.  The concrete
classes are ``Integrator`` and its subclasses ``RK4``, ``RK8``, ``IRKF`` /
``RKF45``, and ``PD45`` (``cpp/lupnt/numerics/integrator.h`` /
``integrator.cc``).

ODE and State Contract
-------------------------------------------------------------------

Every integrator advances a first-order system

.. math::

   \dot{x} = f(t, x),
   \qquad
   x \in \mathbb{R}^n,

where the right-hand side has the LuPNT type

.. math::

   \texttt{ODE} : (t, x) \mapsto f(t, x) \in \mathbb{R}^n .

Numerical dynamics models (e.g. ``NumericalOrbitDynamics``,
``ClockDynamics``, ``JointOrbitClockDynamics``) build an ``ODE`` from their
force/derivative models and pass it to ``Integrator::Propagate`` /
``PropagateEx``.  The integrator does not interpret the state contents; the
unit, frame, and epoch conventions are those of the calling dynamics model
(see :doc:`dynamics`).  The independent variable :math:`t` is in seconds and
the step size :math:`h` is signed by the propagation direction.

Integrator Selection
-------------------------------------------------------------------

``IntegratorType`` selects the concrete subclass instantiated by a dynamics
model's ``SetIntegrator``:

.. math::

   \texttt{IntegratorType} \in
   \{\ \texttt{RK4},\ \texttt{RK8},\ \texttt{RKF45},\ \texttt{PD45}\ \}.

``RK4`` and ``RK8`` are fixed-step methods; ``RKF45`` (an ``IRKF``) and
``PD45`` are embedded, adaptive-step methods.  The default is
``IntegratorType::RK4``.

Fixed-Step Runge-Kutta Contract
-------------------------------------------------------------------

A fixed-step explicit :math:`s`-stage Runge-Kutta method with Butcher
coefficients :math:`(a_{ij}, b_i, c_i)` advances one step of size :math:`h`
by

.. math::

   k_i = f\!\left(t + c_i h,\; x + \sum_{j<i} a_{ij} k_j\right),
   \qquad i = 1,\dots,s,

.. math::

   x_{k+1} = x_k + h \sum_{i=1}^{s} b_i k_i .

In the LuPNT source the stage derivatives are stored pre-scaled by :math:`h`
(``k_i = f(...) * dt``), so the update reads
:math:`x_{k+1} = x_k + \sum_i b_i k_i`.

**RK4** uses the classical tableau
:math:`c = [0,\tfrac12,\tfrac12,1]`, so that

.. math::

   x_{k+1}
   = x_k + \tfrac{1}{6}\left(k_1 + 2k_2 + 2k_3 + k_4\right),

.. math::

   k_1 = h f(t,x),\quad
   k_2 = h f\!\left(t+\tfrac{h}{2}, x+\tfrac{k_1}{2}\right),\quad
   k_3 = h f\!\left(t+\tfrac{h}{2}, x+\tfrac{k_2}{2}\right),\quad
   k_4 = h f(t+h, x+k_3).

**RK8** is a fixed-step 10-stage, 8th-order Runge-Kutta method
(Cooper-Verner-type coefficients).  The stage nodes are
:math:`c = [0,\tfrac{4}{27},\tfrac{2}{9},\tfrac{1}{3},\tfrac12,\tfrac23,
\tfrac16,1,\tfrac56,1]` and the final combination is

.. math::

   x_{k+1}
   = x_k + \tfrac{1}{840}\left(
     41 k_1 + 27 k_4 + 272 k_5 + 27 k_6 + 216 k_7 + 216 k_9 + 41 k_{10}
   \right),

with the intermediate stage coefficients as tabulated in ``RK8::Step``.  Only
the constants that appear in the source are contractual.

Propagation Loop
-------------------------------------------------------------------

``Integrator::Propagate`` marches from :math:`t_0` to :math:`t_f` in steps of
(at most) :math:`\mathrm{d}t`, clamping the final step so it lands exactly on
:math:`t_f`:

.. math::

   h = \min(\mathrm{d}t,\ t_f - t),
   \qquad t \leftarrow t + h .

``PropagateEx`` generalizes this to signed propagation.  With
:math:`\mathrm{dir} = \operatorname{sign}(t_f - t_0)` the step is

.. math::

   h = \mathrm{dir}\cdot\min\!\big(|\mathrm{d}t|,\ |t_f - t|\big),

and :math:`\mathrm{d}t` is treated as a magnitude.  A zero span
(:math:`t_f = t_0`) returns immediately.

Adaptive Step-Size Contract (Embedded RKF)
-------------------------------------------------------------------

``IRKF`` is the abstract embedded Runge-Kutta-Fehlberg base.  Each call to
``Update`` produces a low- and high-order solution pair
:math:`\hat{x}^{\text{low}}_{k+1}`, :math:`\hat{x}^{\text{high}}_{k+1}` for the
same step, whose difference is the local error estimate

.. math::

   e_i = \big|\hat{x}^{\text{high}}_{k+1,i} - \hat{x}^{\text{low}}_{k+1,i}\big| .

Per-component the tolerance is a mixed relative/absolute bound

.. math::

   \mathrm{tol}_i
   = \max\!\big(\texttt{reltol}\cdot |\hat{x}^{\text{high}}_{k+1,i}|,\ \texttt{abstol}\big),
   \qquad
   \varepsilon = \max_i \frac{e_i}{\mathrm{tol}_i}.

The step is **accepted** when :math:`e_i/\mathrm{tol}_i \le
(p+1)^{(p+1)/p}` for every component (where :math:`p` is the low-order method
order; ``ComputeRelError`` uses the non-conservative Butcher acceptance
threshold), and **rejected** otherwise.  The step size is rescaled by a
PI-type controller with safety factor :math:`\beta = 0.9`,

.. math::

   s = \beta \,\varepsilon^{-1/(p+1)},
   \qquad
   s \leftarrow \min\!\big(2.0,\ \max(0.5, s)\big),
   \qquad
   \mathrm{d}t \leftarrow s\,\mathrm{d}t ,

so the step may grow or shrink by at most a factor of two per attempt.
``IRKF::Step`` retries up to ``IntegratorParams::max_iter`` times and returns
the low-order solution; failing to converge throws.

**RKF45** (``IRKF(4)``) is the 6-stage Fehlberg 4(5) pair.  With the stored,
:math:`h`-scaled stages :math:`k_1,\dots,k_6`, the two embedded solutions are

.. math::

   \hat{x}^{\text{low}}_{k+1}
   = x_k + \tfrac{25}{216}k_1 + \tfrac{1408}{2565}k_3
     + \tfrac{2197}{4104}k_4 - \tfrac{1}{5}k_5 ,

.. math::

   \hat{x}^{\text{high}}_{k+1}
   = x_k + \tfrac{16}{135}k_1 + \tfrac{6656}{12825}k_3
     + \tfrac{28561}{56430}k_4 - \tfrac{9}{50}k_5 + \tfrac{2}{55}k_6 ,

with stage nodes :math:`c = [0,\tfrac14,\tfrac38,\tfrac{12}{13},1,\tfrac12]`
and the Fehlberg :math:`a_{ij}` as coded in ``RKF45::Update``.

Adaptive Step-Size Contract (Dormand-Prince)
-------------------------------------------------------------------

``PD45`` is a self-contained Dormand-Prince 4(5) method (it does not use
``IRKF::ComputeRelError``).  It uses the 7-stage tableau ``A_``, the 5th-order
weights ``b_``, and the 4th-order weights ``b_star_`` stored as class
constants.  Per step it forms

.. math::

   y^{\text{high}}_{k+1} = x_k + \sum_{i=1}^{7} b_i\, k_i,
   \qquad
   y^{\text{low}}_{k+1}  = x_k + \sum_{i=1}^{7} b^{*}_i\, k_i ,

with the RMS-scaled error norm

.. math::

   \mathrm{sc}_i = \texttt{abstol} + \texttt{reltol}\,|x_{k,i}|,
   \qquad
   E = \frac{1}{\sqrt{n}}
       \left\|\frac{y^{\text{high}}_{k+1}-y^{\text{low}}_{k+1}}{\mathrm{sc}}\right\|_2 .

The step is accepted when :math:`E \le 1` (returning the high-order solution
:math:`y^{\text{high}}_{k+1}`); otherwise the step is shrunk by

.. math::

   \mathrm{d}t \leftarrow \mathrm{d}t\,\max\!\Big(0.1,\ (0.9/E)^{1/4}\Big),

and retried, up to ``max_iter`` attempts.  The step throws if
:math:`|\mathrm{d}t| < 10^{-3}` or the iteration limit is exceeded.

.. note::

   ``PD45::Step`` evaluates each stage time as :math:`t + a_{i0}\,h` (the
   first column of the tableau) rather than :math:`t + c_i h`; for the
   standard Dormand-Prince nodes these coincide only for the first two
   stages.

Integrator Parameters
-------------------------------------------------------------------

``IntegratorParams`` (set via ``Integrator::SetParams``) carries the adaptive
tolerances and limits, all required positive:

.. math::

   \texttt{max\_iter} \in \mathbb{Z}_{>0},
   \qquad
   \texttt{abstol} > 0,
   \qquad
   \texttt{reltol} > 0 .

It also holds an optional early-termination predicate
:math:`\texttt{terminate\_if}(t, x) \to \{\text{true},\text{false}\}`.

Termination Semantics
-------------------------------------------------------------------

``PropagateEx`` returns an ``IntegratorResult``
:math:`(x, t, \texttt{reason}, \texttt{steps})`.  The predicate is checked
once before stepping and again after every accepted step; the run stops with

.. math::

   \texttt{reason} =
   \begin{cases}
     \texttt{UserCondition} & \text{if } \texttt{terminate\_if}(t,x)=\text{true},\\[2pt]
     \texttt{ReachedTf} & \text{when } t \text{ reaches } t_f .
   \end{cases}

``steps`` counts the accepted integration steps taken.

State-Transition (Sensitivity) Matrix Propagation
-------------------------------------------------------------------

The sensitivity-matrix overloads
``Propagate(..., MatXd* J)`` / ``PropagateEx(..., MatXd* J)`` compute

.. math::

   J = \frac{\partial x(t_f)}{\partial x(t_0)}

by **parallel finite differences** (``JacobianParallel``): the nominal
trajectory is re-propagated for perturbed initial states and the columns of
:math:`J` are assembled from the differences.  This is a finite-difference
sensitivity, not an augmented variational integration; the ordering of rows
and columns of :math:`J` matches the state ordering of :math:`x`.  When ``J``
is ``nullptr`` the Jacobian computation is skipped.

.. note::

   The filter-facing state-transition matrix is obtained separately.
   ``FilterDynamicsFunction`` (see :doc:`filters`) returns the STM
   :math:`F = \partial x(t_f)/\partial x(t_0)` produced by the dynamics
   model's own linearization (autodiff via ``GetFilterDynamicsFunction``),
   which the EKF/SRIF/UDU filters consume directly.  The conceptual
   variational contract :math:`\dot{\Phi} = A(t)\,\Phi`,
   :math:`A = \partial f/\partial x` is documented in :doc:`dynamics`.

Model Boundaries
-------------------------------------------------------------------

* The step size ``dt`` passed to ``Propagate`` must be positive; direction is
  inferred from :math:`t_f - t_0` in ``PropagateEx``.
* ``RK8`` is a fixed-step method with no embedded error estimate; step control
  is available only through ``RKF45`` and ``PD45``.
* The integrator's ``MatXd* J`` sensitivity is finite-difference based and
  therefore re-propagates the trajectory :math:`n+1` times.
