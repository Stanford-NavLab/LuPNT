Measurement Models
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT's non-GNSS
measurement models: inter-satellite crosslinks, ground-station range and
range-rate, lunar surface (LANS/LunaNet) pseudorange, and the terrain-relative
lander sensors (altimeter, crater bearings).  GNSS pseudorange/Doppler/carrier
phase are specified separately in ``docs/math/gnss_measurements.rst``, and the
link budget in ``docs/math/link_budget.rst``.  For each model this page gives
the state dependence, the observation equation :math:`h(x)`, the measurement
Jacobian :math:`H = \partial h/\partial x`, the noise covariance :math:`R`, and
the implementing code.  The concrete classes live under
``cpp/lupnt/measurements/`` and the shared geometry primitives in
``measurement_utils.{h,cc}``.

Measurement Base Contract
-------------------------------------------------------------------

Every state-vector measurement derives from ``Measurement``
(``cpp/lupnt/measurements/measurement.h``) and returns a ``MeasData``
(predicted observable vector :math:`z`, noise covariance :math:`R`) from

.. code-block:: cpp

   virtual MeasData Compute(const State& x, MatXd* H = nullptr) const = 0;

with :math:`z = h(x)`, :math:`R` sized :math:`n_z \times n_z`, and, when
``H`` is non-null, the Jacobian :math:`H = \partial h/\partial x` sized
:math:`n_z \times \dim(x)`.  ``Measurement::CreateFunction`` wraps the model
into a ``FilterMeasurementFunction`` :math:`(x, H, R)\mapsto z` for direct use
with the ``Filter`` family.

Error-state (INS) sensors instead derive from ``ErrorStateMeasurement`` and
consume a ``NavErrorContext`` (nominal position, velocity, attitude
:math:`R_{b\to n}`, clock) rather than a flat ``State``, because their
Jacobian is taken with respect to the nav **error** state
:math:`[\,\delta r,\ \delta v,\ \delta\theta,\ \delta b_a,\ \delta b_g,\
\delta b_\mathrm{clk},\ \delta d_\mathrm{clk}\,]` and depends on the attitude:

.. code-block:: cpp

   virtual MeasData Compute(const NavErrorContext& nominal,
                            MatXd* H = nullptr) const = 0;

Shared Range Primitives
-------------------------------------------------------------------

The geometric building blocks are in
``cpp/lupnt/measurements/measurement_utils.cc``.  For two point masses with
positions :math:`r_1, r_2` and velocities :math:`v_1, v_2` in a common frame,

.. math::

   \rho = \lVert r_2 - r_1 \rVert,
   \qquad
   \dot{\rho}
   =
   \frac{(v_2 - v_1)^\mathsf{T}(r_2 - r_1)}{\lVert r_2 - r_1 \rVert},

.. code-block:: cpp

   Vec2 RangeAndRangeRate(const VecX& r1, const VecX& r2,
                          const VecX& v1, const VecX& v2) {
     Real r_norm = (r2 - r1).norm();
     return Vec2(r_norm, (v2 - v1).dot(r2 - r1) / r_norm);
   }

and the clock-biased forms
:math:`\operatorname{Pseudorange} = \rho + c(b_1 - b_2)`,
:math:`\operatorname{PseudorangeRate} = \dot{\rho} + c(d_1 - d_2)`.

.. note::

   These non-GNSS models use the **instantaneous** geometric range at a single
   epoch: they do not solve the light-time transmit epoch, and they include no
   Shapiro/relativistic/ionospheric corrections (unlike ``GnssMeasurement``).
   The one relativistic term present anywhere in this family is the receiver
   clock scaling :math:`c\,b` in the pseudorange models.

Inter-Satellite Crosslink
-------------------------------------------------------------------

``IslCrosslinkMeasurement``
(``cpp/lupnt/measurements/crosslink_measurement.{h,cc}``) models a hub
satellite's two-way crosslinks to :math:`n_\mathrm{links}` linked satellites,
plus optional one-way pseudoranges to known-position anchor transmitters
(e.g. a lunar surface station with a known clock).

**State dependence.**  The filter state stacks one
:math:`[\,r\ (3),\ v\ (3),\ b_\mathrm{clk},\ d_\mathrm{clk}\,]` block per
satellite (``Config::sub_state_size = 8``): block 0 is the hub, blocks
:math:`1..n_\mathrm{links}` the linked satellites.

**Observation equation.**  Per epoch the model emits, for each link
:math:`i`, a clock-independent two-way range/range-rate pair, followed by one
pseudorange row per anchor :math:`a`:

.. math::

   h(x)
   =
   \begin{bmatrix}
     \big[\,\rho_{0i},\ \dot{\rho}_{0i}\,\big]_{i=1}^{n_\mathrm{links}} \\[4pt]
     \big[\,\lVert r_0 - r_a\rVert + c\,b_0\,\big]_{a}
   \end{bmatrix},

.. code-block:: cpp

   Vec2 yi = RangeAndRangeRate(r0, x.segment(off,3), v0, x.segment(off+3,3));
   y(2*n_links + ad) = (r_hub - ra).norm() + C * b_hub;   // one-way anchor pseudorange

where :math:`\rho_{0i}, \dot{\rho}_{0i}` are the hub-to-satellite-:math:`i`
geometric range/range-rate (``RangeAndRangeRate``), :math:`r_a` an anchor
position (``anchor_pos_mci``, ``Frame::MOON_CI``), and :math:`b_0` the hub
clock bias in seconds.  The two-way crosslink rows are clock-independent
because both directions share the same link.

**Jacobian.**  Taken by autodiff **directly on the raw** :math:`VecX` state:

.. code-block:: cpp

   auto f = [this](const VecX& xx) -> VecX { return Model(xx); };
   jacobian(f, wrt(x_tmp), at(x_tmp), y, *H);

The state is reseeded fresh (dropping incoming derivative seeds), a deliberate
workaround that keeps the seeds intact for the ``SchmidtEKF`` consider-state
update.

**Covariance.**  Diagonal, per row type:

.. math::

   R
   =
   \operatorname{diag}\!\big(
     \underbrace{\sigma_\rho^2, \sigma_{\dot\rho}^2}_{\text{per link}},\ \dots,\
     \underbrace{\sigma_P^2}_{\text{per anchor}}, \dots
   \big),

with defaults :math:`\sigma_\rho = 1` m, :math:`\sigma_{\dot\rho} = 10^{-3}`
m/s, :math:`\sigma_P = 200` m (anchor pseudorange, inflated to absorb the
near-degenerate Earth-GPS-to-surface geometry).

Ground-Station Range and Range-Rate
-------------------------------------------------------------------

``GroundStationRangeMeasurement``
(``cpp/lupnt/measurements/ground_range_measurement.{h,cc}``) models the
instantaneous range and/or range-rate of a target ``[r(3), v(3)]`` relative to
a fixed reference state (the station), with the closed-form design matrix.

**State dependence.**  Target Cartesian state
:math:`x_s = [\,r,\ v\,]`; the station Cartesian state
:math:`x_\mathrm{gs} = [\,r_\mathrm{gs},\ v_\mathrm{gs}\,]` is a fixed
``Config::reference_state`` in the same frame.  With
:math:`\Delta r = r - r_\mathrm{gs}`, :math:`\Delta v = v - v_\mathrm{gs}`,
:math:`u = \Delta r/\rho`.

**Observation equation.**

.. math::

   h(x_s)
   =
   \begin{bmatrix}
     \rho \\
     \dot{\rho}
   \end{bmatrix}
   =
   \begin{bmatrix}
     \lVert \Delta r \rVert \\[2pt]
     \Delta r^\mathsf{T}\Delta v / \rho
   \end{bmatrix}.

Doppler is not emitted as a separate observable; it is the range-rate scaled by
:math:`-f/c` and is obtained from :math:`\dot{\rho}` downstream.

**Jacobian.**  Closed-form (not autodiff), with respect to the target Cart6:

.. code-block:: cpp

   h_obs.block(row,0,1,3) = u.transpose();                        // range
   h_obs.block(row,0,1,3) = ((dv - rho_dot*u)/rho).transpose();   // range-rate, position part
   h_obs.block(row,3,1,3) = u.transpose();                        // range-rate, velocity part

.. math::

   \frac{\partial \rho}{\partial x_s}
   =
   \begin{bmatrix} u^\mathsf{T} & 0 \end{bmatrix},
   \qquad
   \frac{\partial \dot{\rho}}{\partial x_s}
   =
   \begin{bmatrix}
     \dfrac{(\Delta v - \dot{\rho}\,u)^\mathsf{T}}{\rho}
     & u^\mathsf{T}
   \end{bmatrix}.

Batch orbit determination chains this instantaneous design matrix with the
state-transition matrix :math:`\Phi(t_k, t_0)` outside the model; the class
owns only the geometry.

**Covariance.**  Diagonal,
:math:`R = \operatorname{diag}(\sigma_\rho^2, \sigma_{\dot\rho}^2)`, defaults
:math:`\sigma_\rho = 1` m, :math:`\sigma_{\dot\rho} = 10^{-3}` m/s; each row is
optional via ``use_range`` / ``use_range_rate``.

Surface LANS / LunaNet Pseudorange
-------------------------------------------------------------------

``SurfaceLansMeasurement``
(``cpp/lupnt/measurements/surface_measurements.{h,cc}``) is an error-state
model of a one-way code pseudorange from a relay (LCRNS / LANS) satellite to a
lunar surface rover.

**State dependence.**  Error-state
:math:`[\,\delta r,\ \delta v,\ \delta\theta,\ \delta b_a,\ \delta b_g,\
\delta b_\mathrm{clk},\ \delta d_\mathrm{clk}\,]`
(``kSurfaceNavErrorStateSize = 17``); the model reads only the nominal
position :math:`r` and clock bias :math:`b_\mathrm{clk}` from the
``NavErrorContext``.  The transmitter position :math:`r_S` (Moon-fixed) is
supplied with the measurement.

**Observation equation.**

.. math::

   h
   =
   \lVert r_S - r \rVert
   +
   c\,b_\mathrm{clk},

.. code-block:: cpp

   double SurfaceLansMeasurement::PredictedRange(const Vec3d& r_rover,
                                                 double clock_bias_s) const {
     return (r_sat - r_rover).norm() + C * clock_bias_s;
   }

**Jacobian.**  With the rover-to-satellite unit vector
:math:`u = (r_S - r)/\lVert r_S - r\rVert`, the position partial is
:math:`\partial \rho/\partial r = -u`:

.. code-block:: cpp

   H->block(0, config.i_dr, 1, 3) = -u.transpose();   // d(range)/d(r)
   (*H)(0, config.i_dcb) = C;                          // d(rho)/d(clock_bias)

.. math::

   H
   =
   \begin{bmatrix}
     \cdots & -u^\mathsf{T} & \cdots & c & \cdots
   \end{bmatrix},
   \quad
   \frac{\partial h}{\partial \delta r} = -u^\mathsf{T},
   \qquad
   \frac{\partial h}{\partial \delta b_\mathrm{clk}} = c.

**Covariance.**  Scalar; the receiver pseudorange noise and the broadcast
signal-in-space error add in quadrature:

.. math::

   R = \sigma_m^2 + \sigma_\mathrm{sise}^2,

defaults :math:`\sigma_m = 1` m (receiver) and
:math:`\sigma_\mathrm{sise} = 3` m (ephemeris/clock signal-in-space, 1-sigma).
The same model is reused by the lander (``kLanderNavErrorStateSize = 17``).

Lander Altimeter
-------------------------------------------------------------------

``LanderAltimeterMeasurement``
(``cpp/lupnt/measurements/lander_measurements.{h,cc}``) is an error-state model
of a radar/laser altimeter return: the vehicle's height above the terrain
directly below it (nadir), :math:`\text{alt} = U_\mathrm{veh} -
\text{terrain}(E, N)`.

**State dependence.**  Nominal position (error-state offset ``i_dr``).  Because
the predicted altitude and its position partial depend on the DEM, they are
supplied by the simulation (``predicted_altitude_m``, ``h_pos``) and the model
is a thin wrapper:

.. code-block:: cpp

   md.value = VecXd::Constant(1, predicted_altitude_m);
   H->block(0, config.i_dr, 1, 3) = h_pos.transpose();   // d(altitude)/d(r), Moon-fixed

**Observation equation and Jacobian.**

.. math::

   h = \widehat{\text{alt}}(r),
   \qquad
   H
   =
   \begin{bmatrix}
     \cdots & h_\mathrm{pos}^\mathsf{T} & \cdots
   \end{bmatrix},
   \quad
   h_\mathrm{pos}
   =
   \frac{\partial\,\text{alt}}{\partial r}
   \ \text{(DEM slope via local ENU)}.

The DEM-dependent slope couples the vertical measurement to horizontal
position through the local East-North-Up rotation and the terrain gradient
:math:`(\partial\text{terrain}/\partial E,\ \partial\text{terrain}/\partial N)`;
the update is driven by ``LanderNavApp::UpdateAltimeter``.

**Covariance.**  Scalar :math:`R = \sigma_m^2`, default
:math:`\sigma_m = 2` m.

Lander Crater Bearing
-------------------------------------------------------------------

``LanderCraterMeasurement``
(``cpp/lupnt/measurements/lander_measurements.{h,cc}``) is an error-state
terrain-relative-navigation model: the unit line-of-sight from the lander's
down-looking camera to a mapped crater of known Moon-fixed position
:math:`r_c`, expressed in the vehicle **body** frame.

**State dependence.**  Nominal position :math:`r` and attitude
:math:`R_{b\to n}`; error-state offsets ``i_dr`` (position) and ``i_dth``
(attitude).  The residual constrains both horizontal position and attitude.

**Observation equation.**  With the nav-frame line-of-sight
:math:`u_n = (r_c - r)/\lVert r_c - r\rVert`,

.. math::

   h = R_{b\to n}^\mathsf{T}\, u_n,
   \qquad
   u_n = \frac{r_c - r}{\lVert r_c - r\rVert}.

.. code-block:: cpp

   Vec3d u_n = los / range;                 // nav-frame unit LOS
   md.value = nom.R_b2n.transpose() * u_n;  // predicted body LOS

**Jacobian.**  Coupling position and attitude, with
:math:`\partial u_n/\partial r = -\tfrac{1}{\rho}(I - u_n u_n^\mathsf{T})`:

.. code-block:: cpp

   Mat3d dun_dr = -(1.0/range) * (Mat3d::Identity() - u_n * u_n.transpose());
   H->block(0, config.i_dr,  3, 3) = nom.R_b2n.transpose() * dun_dr;
   H->block(0, config.i_dth, 3, 3) = nom.R_b2n.transpose() * Skew3d(u_n);

.. math::

   \frac{\partial h}{\partial \delta r}
   =
   R_{b\to n}^\mathsf{T}
   \left[-\frac{1}{\rho}\big(I - u_n u_n^\mathsf{T}\big)\right],
   \qquad
   \frac{\partial h}{\partial \delta\theta}
   =
   R_{b\to n}^\mathsf{T}\,[u_n]_\times,

where :math:`[u_n]_\times` is the skew-symmetric cross-product matrix.

**Covariance.**  Isotropic per line-of-sight component (small-angle
unit-vector approximation),

.. math::

   R = \sigma_\mathrm{rad}^2\, I_3,

default :math:`\sigma_\mathrm{rad} = 10^{-3}` rad.  A degenerate zero-range
geometry returns an empty value vector, so the hosting application skips the
update.

Model Boundaries
-------------------------------------------------------------------

* All non-GNSS models use instantaneous geometric range; no light-time
  iteration and no Shapiro/relativistic/ionospheric range corrections (the
  only relativistic term is the :math:`c\,b` clock scaling in the pseudorange
  models).
* ``IslCrosslinkMeasurement`` reseeds a fresh autodiff state to preserve
  consider-state seeds for ``SchmidtEKF``; ``GroundStationRangeMeasurement``
  uses a closed-form design matrix.
* Error-state models (``SurfaceLansMeasurement``, ``LanderAltimeter``,
  ``LanderCrater``) linearize about a ``NavErrorContext`` and drive the
  Joseph-form error-state update from the hosting application rather than
  through ``FilterMeasurementFunction``.
* The altimeter's predicted value and Jacobian are DEM-supplied; the model
  itself carries no terrain evaluation.
* Covariances are diagonal (or isotropic) by convention.
