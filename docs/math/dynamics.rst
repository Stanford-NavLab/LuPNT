Dynamics Models
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT orbit dynamics,
with emphasis on the numerical Cartesian dynamics used by
``NBodyDynamics``.  The goal is to make the state, epoch, unit, frame, and
force-model conventions explicit enough that propagation and filtering tests
can compare against the same equations.

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

Setting the body-fixed frame equal to the integration frame recovers the older
inertial-axis J2 convention and is useful for reference comparisons that were
generated with that approximation.

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

Relativistic Orbit Correction
-------------------------------------------------------------------

When relativity is enabled, ``NBodyDynamics`` applies this correction for the
Sun and for the closest non-Sun, non-SSB configured body, interpreted as the
local central planet.  For the central body, let

.. math::

   r_c = r - r_B,
   \qquad
   v_c = v - v_B,
   \qquad
   \mu = \mu_B.

The first-order Schwarzschild post-Newtonian acceleration correction is

.. math::

   a_\mathrm{rel}
   =
   -\frac{\mu}{c^2\lVert r_c\rVert^3}
   \left[
     \left(
       \frac{4\mu}{\lVert r_c\rVert}
       -
       v_c^\mathsf{T}v_c
     \right)r_c
     +
     4(r_c^\mathsf{T}v_c)v_c
   \right].

This is the gravito-electric correction.  Frame-dragging and third-body
post-Newtonian terms are not included in the current dynamics model.

The Sun term uses the same expression with

.. math::

   r_c = r - r_\odot,
   \qquad
   v_c = v - v_\odot,
   \qquad
   \mu = \mu_\odot.

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
