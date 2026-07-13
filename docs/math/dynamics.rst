Dynamics Models
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT orbit dynamics,
with emphasis on the numerical Cartesian dynamics used by
``NBodyDynamics``.  The goal is to make the state, epoch, unit, frame, and
force-model conventions explicit enough that propagation and filtering tests
can compare against the same equations.

The Cartesian rate assembly lives in
``cpp/lupnt/dynamics/numerical_orbit_dynamics.cc`` (classes
``NBodyDynamics``, ``JToCartTwoBodyDynamics``, ``CartesianTwoBodyDynamics``);
the individual force-model accelerations live in
``cpp/lupnt/environment/forces.cc`` / ``forces.h``.

Unified Force-Model Configuration
-------------------------------------------------------------------

Every scenario specifies its orbit force model with the same YAML ``force_model:``
block, whether it is the shared truth model (``world.force_model``), an agent's
``dynamics:``, or an estimator's filter/truth dynamics.  The block is a ``bodies:``
list plus optional ``relativity`` and cannonball-SRP (``CR`` / ``area`` / ``mass``)
entries:

.. code-block:: yaml

   force_model:
     bodies:
       - MOON: {n: 20, m: 20}   # spherical-harmonic gravity field to degree/order 20
       - EARTH: {}              # third body, point-mass (empty {} = no harmonics)
       - SUN: {}                # third body, point-mass
     relativity: true           # n-body post-Newtonian correction (see below)
     CR: 1.0                    # SRP: enables cannonball SRP with B_SRP = CR * area / mass
     area: 0.002
     mass: 1.0

Each ``bodies`` entry is ``BODY: {n: <degree>, m: <order>}`` for a
spherical-harmonic gravity field; an empty ``{}`` (or omitted ``n``/``m``) selects
point-mass gravity.  A body's presence in the list is what enables its
contribution to the sum in `Total N-Body Acceleration`_ below.

The same block feeds two consumers so that world, agent, and filter dynamics stay
consistent:

* ``NBodyDynamics(Config&)``
  (``cpp/lupnt/dynamics/numerical_orbit_dynamics.cc``) constructs a full dynamics
  object directly from the block — adding each body via ``Body::CreateBody``,
  setting the SRP/drag ballistic coefficients from ``CR``/``CD``, ``area``,
  ``mass``, and reading ``frame``, ``units``, ``autodiff``, and ``relativity``.
* ``ParseForceModelSpec(const Config&)`` reads the same block into a lightweight
  ``ForceModelSpec`` (``moon_degree``/``moon_order`` from the ``MOON`` entry,
  ``include_earth``/``include_sun`` from the presence of those bodies, plus
  ``relativity`` and the SRP scalars) so per-example dynamics builders that need
  their own integrator tolerances or clock coupling can be fed from one
  consistent config surface.  Both ``relativity`` and ``use_relativity`` are
  accepted.

State, Epoch, Frame, and Unit Contract
-------------------------------------------------------------------

The propagated Cartesian orbit state is

.. math::

   x =
   \begin{bmatrix}
     r \\
     v
   \end{bmatrix},
   \qquad
   r \in \mathbb{R}^3,
   \qquad
   v \in \mathbb{R}^3.

All components of :math:`r`, :math:`v`, and :math:`\dot{x}` are expressed in
the configured integration frame ``frame_`` and configured coherent
``UnitSystem``:

.. math::

   [r] = L, \qquad
   [v] = L/T, \qquad
   [a] = L/T^2.

For ``NBodyDynamics``, the integration variable :math:`t` is a simulation time
offset.  The ephemeris epoch used internally is

.. math::

   t_\mathrm{TDB} = t_\mathrm{LuPNT,0} + t,

where :math:`t_\mathrm{LuPNT,0}` is ``GetLupntEpoch()``.  Ephemeris and frame
queries are therefore evaluated at TDB seconds from the LuPNT epoch origin.

The model computes

.. math::

   \dot{x}
   =
   f(t, x)
   =
   \begin{bmatrix}
     v \\
     a(t, r, v)
   \end{bmatrix}.

Implemented by
``cpp/lupnt/dynamics/numerical_orbit_dynamics.cc :: NBodyDynamics::ComputeRates``.
The TDB epoch offset and the :math:`\dot{x} = [v;\,a]` assembly are the first
and last lines of the routine:

.. code-block:: cpp

   VecX NBodyDynamics::ComputeRates(Real t, const State& rv) const {
     Real t_tdb = t + GetLupntEpoch();
     Vec3 r = rv.head(3);
     Vec3 v = rv.tail(3);
     Vec3 a = Vec3::Zero();
     // ... sum force-model accelerations into a ...
     Vec6 rv_dot;
     rv_dot << v, a;
     return rv_dot;
   }

Unit Conversion
-------------------------------------------------------------------

Historical LuPNT constants are SI-valued.  When a model uses a non-SI
``UnitSystem``, dimensional values are converted through

.. math::

   q_\mathrm{unit}
   =
   \frac{q_\mathrm{SI}}
        {L^p T^q M^s},

where :math:`p`, :math:`q`, and :math:`s` are the length, time, and mass
powers of the quantity.

For example,

.. math::

   r_\mathrm{unit} = \frac{r_\mathrm{SI}}{L},
   \qquad
   v_\mathrm{unit} = \frac{v_\mathrm{SI}}{L/T},
   \qquad
   \mu_\mathrm{unit} = \frac{\mu_\mathrm{SI}}{L^3/T^2}.

Planetary ephemerides are evaluated in the internal TDB-compatible coordinate
scale and then converted to the requested unit system at the dynamics boundary.

Implemented by the unit-conversion helpers in the anonymous namespace of
``cpp/lupnt/dynamics/numerical_orbit_dynamics.cc``
(``PositionFromSI`` / ``PositionToSI`` / ``VelocityToSI`` /
``AccelerationFromSI`` / ``StateToSI``), each a scaling by the ``UnitSystem``
length/time powers:

.. code-block:: cpp

   Vec3 PositionFromSI(const Vec3& r, const UnitSystem& units) { return r / units.length; }
   Vec3 VelocityToSI(const Vec3& v, const UnitSystem& units) {
     return v * (units.length / units.time);
   }
   Vec3 AccelerationFromSI(const Vec3& a, const UnitSystem& units) {
     return a * (units.time * units.time / units.length);
   }

Point-Mass Gravity
-------------------------------------------------------------------

For a perturbing body with gravitational parameter :math:`\mu_i` and position
:math:`s_i` in the integration frame, LuPNT uses the relative point-mass
acceleration

.. math::

   a_i
   =
   -\mu_i
   \left(
     \frac{r - s_i}{\lVert r - s_i \rVert^3}
     +
     \frac{s_i}{\lVert s_i \rVert^3}
   \right).

The second term removes the acceleration of the integration-frame origin when
the frame origin is not the perturbing body.  If :math:`s_i = 0`, only the
central two-body term remains:

.. math::

   a_i
   =
   -\mu_i \frac{r}{\lVert r \rVert^3}.

Implemented by
``cpp/lupnt/environment/forces.cc :: AccelerationPointMass`` (called per
non-gravity-field body by ``NBodyDynamics::ComputeRates``):

.. code-block:: cpp

   Vec3 AccelerationPointMass(const Vec3& r, const Vec3& s, Real GM) {
     Vec3 d = r - s;
     Vec3 a = Vec3::Zero();
     if (s.norm() > EPS) a += s / pow(s.norm(), 3);   // indirect (origin recoil) term
     if (d.norm() > EPS) a += d / pow(d.norm(), 3);   // direct attraction
     a *= -GM;
     return a;
   }

Gravity-Field Acceleration
-------------------------------------------------------------------

For a body with a configured spherical-harmonic gravity field, the spacecraft
position is first rotated into the body-fixed frame:

.. math::

   r_\mathrm{bf}
   =
   R_{\mathrm{frame}\rightarrow\mathrm{bf}}(t_\mathrm{TDB})\,r.

The gravity field acceleration is computed in that fixed frame from
unnormalized coefficients :math:`C_{nm}` and :math:`S_{nm}` up to configured
degree and order:

.. math::

   a_\mathrm{bf}
   =
   \nabla U(r_\mathrm{bf}; \mu, R_\mathrm{ref}, C_{nm}, S_{nm}).

The result is rotated back to the integration frame:

.. math::

   a_\mathrm{field}
   =
   R_{\mathrm{bf}\rightarrow\mathrm{frame}}(t_\mathrm{TDB})\,a_\mathrm{bf}.

This path uses only the rotation component of the frame transform for
accelerations.

Implemented by
``cpp/lupnt/environment/forces.cc :: AccelarationGravityField`` (the
Montenbruck-Gill ``V_nm`` / ``W_nm`` harmonic recursion, templated on ``Real``
for autodiff or ``double`` for fast propagation), driven by the gravity-field
branch of ``NBodyDynamics::ComputeRates``:

.. code-block:: cpp

   // NBodyDynamics::ComputeRates -- rotate to body-fixed, evaluate field, rotate back
   Vec3 r_bf = PositionFromSI(ConvertFrame(t_tdb, r_si, frame_, body.fixed_frame), units_);
   a_bf = AccelarationGravityField<Real>(r_bf, grav.GM, grav.R, grav.CS, grav.n, grav.m);
   auto [R_bf_to_frame, translation] = GetFrameRotationTranslation(t_tdb, body.fixed_frame, frame_);
   a += R_bf_to_frame * a_bf;

J2 Cartesian Dynamics
-------------------------------------------------------------------

The Cartesian J2 model augments two-body gravity with a zonal term evaluated
about the body's fixed spin axis, not about the integration-frame z-axis.

Let :math:`F` be the integration frame and :math:`B` be the body-fixed frame
configured on ``JToCartTwoBodyDynamics``.  These frames must share the same
origin.  The spacecraft position is first rotated into the body-fixed frame:

.. math::

   r_B = R_{BF}(t_\mathrm{TDB}) r_F.

The J2 acceleration is computed in the body-fixed frame:

.. math::

   a_{J2,x}
   =
   -\frac{3}{2}
   \frac{\mu J_2 R^2}{r^5}
   \left(1 - 5\frac{z^2}{r^2}\right)x,

.. math::

   a_{J2,y}
   =
   -\frac{3}{2}
   \frac{\mu J_2 R^2}{r^5}
   \left(1 - 5\frac{z^2}{r^2}\right)y,

.. math::

   a_{J2,z}
   =
   -\frac{3}{2}
   \frac{\mu J_2 R^2}{r^5}
   \left(3 - 5\frac{z^2}{r^2}\right)z,

where :math:`[x,y,z]^T = r_B` and
:math:`r = \lVert r_B \rVert = \lVert r_F \rVert`.

The acceleration is then rotated back to the integration frame:

.. math::

   a_{J2,F}
   =
   R_{BF}^\mathsf{T} a_{J2,B}.

Setting the body-fixed frame equal to the integration frame yields the
inertial-axis J2 convention, useful for reference comparisons that use that
approximation.

Implemented by
``cpp/lupnt/dynamics/numerical_orbit_dynamics.cc :: JToCartTwoBodyDynamics::ComputeRates``:

.. code-block:: cpp

   Vec3 r_bf = R_frame_to_body_fixed * r;
   Real aux1 = -3.0 / 2.0 * GM_ * J2_ * pow(R_body_, 2.0) / pow(r_norm, 5.0);
   Real aux2 = 5.0 * pow(r_bf(2) / r_norm, 2.0);
   a_J2_bf(0) = aux1 * (1.0 - aux2) * r_bf(0);
   a_J2_bf(1) = aux1 * (1.0 - aux2) * r_bf(1);
   a_J2_bf(2) = aux1 * (3.0 - aux2) * r_bf(2);
   rv_dot.tail(3) += R_frame_to_body_fixed.transpose() * a_J2_bf;

The plain two-body term :math:`-\mu r/\lVert r\rVert^3` is
``cpp/lupnt/dynamics/numerical_orbit_dynamics.cc :: CartesianTwoBodyDynamics::ComputeRates``.

Solar Radiation Pressure
-------------------------------------------------------------------

Solar radiation pressure is enabled by ``SetUseSrp`` or by setting the SRP
ballistic coefficient

.. math::

   B_\mathrm{SRP} = C_R \frac{A}{m}.

For a Sun vector :math:`r_\odot` in the integration frame, the SRP acceleration
is

.. math::

   a_\mathrm{SRP}
   =
   \nu(r, r_\odot)
   B_\mathrm{SRP}
   P_\odot
   \mathrm{AU}^2
   \frac{r - r_\odot}{\lVert r - r_\odot \rVert^3}.

The illumination factor :math:`\nu` is the apparent-disk overlap shadow
function:

.. math::

   0 \le \nu \le 1,

with :math:`\nu = 0` in umbra, :math:`\nu = 1` in full sunlight, and partial
values in penumbra.

Let

.. math::

   \alpha = \sin^{-1}\left(\frac{R_\odot}{\lVert r_\odot-r\rVert}\right),
   \qquad
   \beta = \sin^{-1}\left(\frac{R_B}{\lVert r\rVert}\right),

.. math::

   \gamma =
   \cos^{-1}
   \left(
     \frac{-r^\mathsf{T}(r_\odot-r)}
          {\lVert r\rVert \lVert r_\odot-r\rVert}
   \right).

The spacecraft is fully illuminated when

.. math::

   \gamma \ge \alpha + \beta.

It is in umbra when the occulting disk covers the solar disk:

.. math::

   \gamma \le |\beta-\alpha|,
   \qquad
   \beta \ge \alpha.

Otherwise, :math:`\nu` is one minus the overlap area of the two apparent disks
divided by the solar disk area.

The cannonball acceleration is
``cpp/lupnt/environment/forces.cc :: AccelerationSolarRadiation`` and the
shadow factor :math:`\nu` is
``cpp/lupnt/environment/forces.cc :: ShadowFunction`` (wrapped by
``Illumination``); ``NBodyDynamics::ComputeRates`` multiplies them together.
The two apparent-radius angles and the umbra / penumbra tests are:

.. code-block:: cpp

   // ShadowFunction: apparent radii of Sun (a) and occulting body (b), separation (c)
   Real a = asin(ClampUnit(R_sun / rho_sun_norm));
   Real b = asin(ClampUnit(R_body / r_norm));
   Real c = acos(ClampUnit(-r.dot(rho_sun) / (r_norm * rho_sun_norm)));
   if (c >= a + b) return 1.0;                       // fully illuminated
   if (c <= abs(b - a)) { if (b >= a) return 0.0; }  // umbra / annular

.. note::

   ``ClampUnit`` returns a strictly-interior constant :math:`\pm(1-10^{-12})`
   at the domain edges of ``asin`` / ``acos`` so the autodiff derivative stays
   finite across a grazing shadow boundary; this is what lets the SRP
   state-transition matrix be formed analytically (see the ``ClampUnit`` /
   ``Clamp01`` comments in ``forces.cc``).

Atmospheric Drag
-------------------------------------------------------------------

Earth atmospheric drag uses a Harris-Priester density model.  The ballistic
coefficient is

.. math::

   B_D = C_D \frac{A}{m}.

In the true-of-date frame,

.. math::

   v_\mathrm{rel}
   =
   v_\mathrm{tod} - \omega_E \times r_\mathrm{tod},

and the drag acceleration is

.. math::

   a_D
   =
   -\frac{1}{2}
   B_D
   \rho_\mathrm{HP}
   \lVert v_\mathrm{rel} \rVert
   v_\mathrm{rel}.

The implementation converts the propagated state to SI for the drag model and
then converts the resulting acceleration back to the configured unit system.

Implemented by
``cpp/lupnt/environment/forces.cc :: AccelerationDrag`` (with the density from
``cpp/lupnt/environment/forces.cc :: DensityHarrisPriester``), invoked for
Earth by ``NBodyDynamics::ComputeRates``:

.. code-block:: cpp

   Vec3 v_rel = v_tod - omega.cross(r_tod);          // Earth-relative velocity
   Real v_abs = v_rel.norm();
   Real dens  = DensityHarrisPriester(mjd_tt, r_tod);
   Vec3 a_tod = -0.5 * bcoeff_drag * dens * v_abs * v_rel * KM_M;
   return T.transpose() * a_tod;                      // back to inertial

``DensityHarrisPriester`` interpolates the tabulated min/max density profiles
exponentially in altitude and weights them by the diurnal-bulge factor
:math:`\cos^{n}(\psi/2)` before returning kg/m^3.

Relativistic Orbit Correction
-------------------------------------------------------------------

When ``relativity`` is enabled, ``NBodyDynamics`` adds the full *n*-body
point-mass relativistic perturbative acceleration of Moyer (2000), Eq. (4-26) --
the parameterized post-Newtonian (PPN) Einstein--Infeld--Hoffmann acceleration in
the Solar-System barycentric frame.  The leading Newtonian point-mass term (the
``1`` in Moyer's first brace) is removed so this quantity *adds* to the Newtonian
gravity already summed above.

Let :math:`i` denote the spacecraft and :math:`j,k,l` the configured massive
bodies, with barycentric positions :math:`r`, velocities :math:`\dot r`,
gravitational parameters :math:`\mu`, and pairwise distances
:math:`r_{ij}=\lVert r_i-r_j\rVert`.  With PPN parameters :math:`\beta` and
:math:`\gamma` (both unity in general relativity) and speed of light :math:`c`,

.. math::

   \begin{aligned}
   a_\mathrm{rel}
   ={}& \sum_{j\neq i}\frac{\mu_j\,(r_j-r_i)}{r_{ij}^3}
     \Big\{
       -\frac{2(\beta+\gamma)}{c^2}\sum_{l\neq i}\frac{\mu_l}{r_{il}}
       -\frac{2\beta-1}{c^2}\sum_{k\neq j}\frac{\mu_k}{r_{jk}} \\
   &\qquad\quad
       +\gamma\frac{\lVert\dot r_i\rVert^2}{c^2}
       +(1+\gamma)\frac{\lVert\dot r_j\rVert^2}{c^2}
       -\frac{2(1+\gamma)}{c^2}\,\dot r_i\!\cdot\!\dot r_j \\
   &\qquad\quad
       -\frac{3}{2c^2}\!\left[\frac{(r_i-r_j)\!\cdot\!\dot r_j}{r_{ij}}\right]^2
       +\frac{1}{2c^2}(r_j-r_i)\!\cdot\!\ddot r_j
     \Big\} \\
   &+\frac{1}{c^2}\sum_{j\neq i}\frac{\mu_j}{r_{ij}^3}
     \Big\{(r_i-r_j)\!\cdot\!\big[(2+2\gamma)\dot r_i-(1+2\gamma)\dot r_j\big]\Big\}
     (\dot r_i-\dot r_j) \\
   &+\frac{3+4\gamma}{2c^2}\sum_{j\neq i}\frac{\mu_j\,\ddot r_j}{r_{ij}}.
   \end{aligned}

The perturbing-body accelerations :math:`\ddot r_j` are taken from the Newtonian
*n*-body model, :math:`\ddot r_j=\sum_{k\neq j}\mu_k(r_k-r_j)/r_{jk}^3`; terms of
order :math:`1/c^4` are dropped, so this Newtonian value is sufficient (Moyer,
p. 4-21).

Because the equation is written in barycentric coordinates, position
*differences* are frame-independent but the absolute velocities in the
:math:`1/c^2` terms are not.  ``NBodyDynamics`` therefore gathers
Solar-System-barycenter (``BodyId::SSB``) referenced states of the spacecraft and
of every configured body before evaluating it.

For a single perturbing body at rest, the model above reduces *exactly* to the
one-body Schwarzschild isotropic form (Moyer Eq. (4-61)),

.. math::

   a_\mathrm{rel}
   =
   \frac{\mu}{c^2\lVert r\rVert^3}
   \left[
     \left(2(\beta+\gamma)\frac{\mu}{\lVert r\rVert}
           -\gamma\lVert\dot r\rVert^2\right)r
     +2(1+\gamma)(r\!\cdot\!\dot r)\,\dot r
   \right],

which is the regression check in
``cpp/test/dynamics/test_relativity_nbody.cc``.  Frame-dragging (Lense--Thirring)
and geodesic-precession terms are not included in the current dynamics model.

Implemented by
``cpp/lupnt/environment/forces.cc :: AccelerationRelativisticNBody`` and wired
into ``NBodyDynamics::ComputeRates`` and ``ComputeAccelerations`` through the
private ``RelativisticNBodyAcceleration`` helper, which assembles the
SSB-referenced body states.

.. note::

   Reference: T. D. Moyer, *Formulation for Observed and Computed Values of Deep
   Space Network Data Types for Navigation*, JPL Publication 00-7 (DESCANSO
   Monograph 2), 2000, Eqs. (4-26) and (4-61).

Total N-Body Acceleration
-------------------------------------------------------------------

The total acceleration assembled by ``NBodyDynamics`` is

.. math::

   a
   =
   \sum_i a_i
   +
   \sum_j a_{\mathrm{field},j}
   +
   a_\mathrm{SRP}
   +
   a_D
   +
   a_\mathrm{rel},

with optional terms omitted when their corresponding switches are disabled or
when no applicable body is configured.

Assembled by
``cpp/lupnt/dynamics/numerical_orbit_dynamics.cc :: NBodyDynamics::ComputeRates``
(the per-body loop adds gravity, SRP, and drag; the relativity terms are added
after the loop).  For a per-term breakdown of the same sum, see
``cpp/lupnt/dynamics/numerical_orbit_dynamics.cc :: NBodyDynamics::ComputeAccelerations``,
which returns a map keyed ``"<BODY>_gravity"``, ``"srp"``, ``"drag"``,
``"relativity"``, etc., whose values sum to the acceleration part of
``ComputeRates``.

State Transition Matrix Convention
-------------------------------------------------------------------

When STM propagation is requested, the integrator propagates the variational
system for the same state ordering:

.. math::

   \Phi(t_f,t_0)
   =
   \frac{\partial x(t_f)}{\partial x(t_0)}.

For a Cartesian state,

.. math::

   \dot{\Phi}
   =
   A(t)\Phi,
   \qquad
   A(t)
   =
   \frac{\partial f}{\partial x}.

The finite-dimensional contract is that rows and columns follow the order

.. math::

   [r_x,\ r_y,\ r_z,\ v_x,\ v_y,\ v_z].

Requested through
``cpp/lupnt/dynamics/numerical_orbit_dynamics.cc :: NBodyDynamics::Propagate``
(the ``MatXd* stm`` overload), which requires autodiff to be enabled:

.. code-block:: cpp

   State NBodyDynamics::Propagate(const State& x0, Real t0, Real tf, const State* u, MatXd* stm) {
     LUPNT_CHECK(use_ad_, "Autodiff not enabled", "NBodyDynamics");
     return NumericalDynamics::Propagate(x0, t0, tf, u, stm);
   }

See :doc:`integration` for how the sensitivity matrix is actually formed
(``Integrator::Propagate(..., MatXd* J)`` via ``JacobianParallel``).

Current Model Boundaries
-------------------------------------------------------------------

The following are intentional boundaries of the current implementation:

* ``NBodyDynamics`` expects its configured bodies to use the same
  ``UnitSystem`` as the dynamics model.
* ``SetUnits`` must be called before adding bodies.
* Relativity is evaluated from the closest configured planet-like body, not
  from every gravitating body.
* SRP uses a cannonball coefficient and spherical occulting bodies.
* Drag is currently Earth-specific.
* Force-model accelerations are expressed in the integration frame and
  configured unit system before being returned by ``ComputeRates``.
