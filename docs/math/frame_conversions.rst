Frame Conversions
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT frame
conversion functions such as ``ConvertFrame``,
``GetFrameRotationTranslation``, and the Earth/Moon frame-specific conversion
helpers.  The key requirements are explicit epoch scale, frame origin,
rotation, translation, and velocity conventions.

The public dispatch (``ConvertFrame``, ``GetFrameRotationTranslation``,
``GetFrameRotationTranslationRv``, ``GetFrameCenter``) lives in
``cpp/lupnt/conversions/frame_converter.cc`` /
``frame_converter.h``; the individual body-pair rotations
(``GcrfToItrf``, ``MoonCiToPa``, the IAU planet orientation, ...) live in
``cpp/lupnt/conversions/frame_conversions.cc`` / ``frame_conversions.h``.

Epoch Contract
-------------------------------------------------------------------

Frame conversions use epochs expressed as TDB seconds from the LuPNT epoch
origin:

.. math::

   t \equiv t_\mathrm{TDB}.

If an input time is in another time scale, it must be converted to TDB before
calling frame conversion functions.  Earth-orientation internals convert from
TDB to TT, UT1, or UTC as needed to evaluate precession-nutation, Earth
rotation, polar motion, and EOP corrections.

Affine Position Transform
-------------------------------------------------------------------

For position-only conversions, LuPNT represents the mapping from frame
:math:`A` to frame :math:`B` as an affine transform:

.. math::

   r_B
   =
   R_{BA}(t) r_A
   +
   p_{BA}(t).

Here :math:`R_{BA}` rotates vector components from frame :math:`A` axes into
frame :math:`B` axes, and :math:`p_{BA}` is the position of the origin of
frame :math:`A` expressed in frame :math:`B`.

``GetFrameRotationTranslation`` returns exactly the pair
:math:`(R_{BA}, p_{BA})`.

Implemented by
``cpp/lupnt/conversions/frame_converter.cc :: GetFrameRotationTranslation``,
which recovers :math:`(R, p)` by transforming the identity basis and the zero
vector through ``ConvertFrame``:

.. code-block:: cpp

   Mat3 M = ConvertFrame(t_tdb, I3, from_frame, to_frame);          // basis images + offset
   Vec3 r = ConvertFrame(t_tdb, MatX3::Zero(1, 3), from_frame, to_frame).transpose();
   Mat3 R = (M.rowwise() - r.transpose()).transpose();              // subtract offset -> rotation
   return std::make_pair(R, r);

Cartesian State Transform
-------------------------------------------------------------------

For a Cartesian state, the full transform is

.. math::

   \begin{bmatrix}
     r_B \\
     v_B
   \end{bmatrix}
   =
   R_{6,BA}(t)
   \begin{bmatrix}
     r_A \\
     v_A
   \end{bmatrix}
   +
   p_{6,BA}(t).

For a pure time-dependent rotation and translation,

.. math::

   r_B
   =
   R_{BA} r_A + p_{BA},

.. math::

   v_B
   =
   R_{BA} v_A
   +
   \dot{R}_{BA} r_A
   +
   \dot{p}_{BA}.

``GetFrameRotationTranslationRv`` returns the affine six-state equivalent
:math:`(R_{6,BA}, p_{6,BA})` by applying ``ConvertFrame`` to the Cartesian
basis states and the zero state.

Implemented by
``cpp/lupnt/conversions/frame_converter.cc :: GetFrameRotationTranslationRv``:

.. code-block:: cpp

   Vec6 t6 = ConvertFrame(t_tdb, MatX6::Zero(1, 6), from_frame, to_frame, false).row(0).transpose();
   MatX6 Y = ConvertFrame(t_tdb, Mat6::Identity(), from_frame, to_frame, false);   // 6x6
   for (int j = 0; j < 6; ++j) R6.col(j) = Y.row(j).transpose() - t6;              // columns of R6
   return {R6, t6};

Rotate-Only Mode
-------------------------------------------------------------------

When ``rotate_only`` is true, ``ConvertFrame`` applies only the rotation
returned by ``GetFrameRotationTranslation``:

.. math::

   r_B = R_{BA}r_A + p_{BA},
   \qquad
   v_B = R_{BA}v_A.

This mode is useful for vector-like quantities such as accelerations or local
offsets where the origin velocity terms should not be applied.  For ordinary
position-velocity state transformations, ``rotate_only`` should be false.

Implemented by the ``rotate_only`` branch of the ``Vec6`` overload of
``cpp/lupnt/conversions/frame_converter.cc :: ConvertFrame``:

.. code-block:: cpp

   Vec6 ConvertFrame(Real t_tdb, const Vec6& rv_in, Frame frame_in, Frame frame_out, bool rotate_only) {
     if (!rotate_only) return ConvertFrameBase(t_tdb, rv_in, frame_in, frame_out);
     auto [R, r] = GetFrameRotationTranslation(t_tdb, frame_in, frame_out);
     Vec6 rv_out;
     rv_out.head(3) = R * rv_in.head(3) + r;   // position: rotate + translate
     rv_out.tail(3) = R * rv_in.tail(3);       // velocity: rotate only
     return rv_out;
   }

Frame Graph
-------------------------------------------------------------------

The implemented native conversion graph is centered on GCRF and MOON_CI:

.. math::

   \mathrm{ICRF}
   \leftrightarrow
   \mathrm{GCRF}
   \leftrightarrow
   \mathrm{MOON\_CI}.

Earth-fixed and Earth inertial frames attach through GCRF:

.. math::

   \mathrm{ITRF}
   \leftrightarrow
   \mathrm{GCRF}
   \leftrightarrow
   \mathrm{EME}.

Lunar fixed and lunar orbit-plane frames attach through MOON_CI:

.. math::

   \mathrm{MOON\_ME}
   \leftrightarrow
   \mathrm{MOON\_PA}
   \leftrightarrow
   \mathrm{MOON\_CI}
   \leftrightarrow
   \mathrm{MOON\_OP}.

The Earth-Moon rotating frame EMR attaches through GCRF.

The routing itself is
``cpp/lupnt/conversions/frame_converter.cc :: ConvertFrameBase`` (a recursive
``switch`` on ``frame_in`` that hops through the GCRF / MOON_CI / ICRF hubs);
solar-system planet frames are split off first to
``ConvertPlanetFrame``:

.. code-block:: cpp

   // ConvertFrameBase -- e.g. GCRF -> any lunar frame routes through MOON_CI
   case Frame::GCRF:
     switch (frame_out) {
       case Frame::ITRF: return GcrfToItrf(t_tdb, rv_in);
       case Frame::MOON_PA: case Frame::MOON_ME: case Frame::MOON_OP: case Frame::MOON_CI:
         return ConvertFrame(t_tdb, GcrfToMoonCi(t_tdb, rv_in), Frame::MOON_CI, frame_out);
     }

Frame Centers
-------------------------------------------------------------------

Each supported frame has a natural center used by ephemeris and conversion
helpers:

.. math::

   \operatorname{center}(\mathrm{ICRF}) = \mathrm{SSB},

.. math::

   \operatorname{center}(\mathrm{GCRF})
   =
   \operatorname{center}(\mathrm{ITRF})
   =
   \operatorname{center}(\mathrm{EME})
   =
   \mathrm{Earth},

.. math::

   \operatorname{center}(\mathrm{MOON\_CI})
   =
   \operatorname{center}(\mathrm{MOON\_PA})
   =
   \operatorname{center}(\mathrm{MOON\_ME})
   =
   \operatorname{center}(\mathrm{MOON\_OP})
   =
   \mathrm{Moon}.

Implemented by the ``frame_centers`` map and
``cpp/lupnt/conversions/frame_converter.cc :: GetFrameCenter`` (also covers the
solar-system ``<PLANET>_CI`` / ``<PLANET>_FIXED`` frames).

GCRF and ICRF
-------------------------------------------------------------------

GCRF and ICRF are related by the Earth barycentric state.  If
:math:`r_{\mathrm{SSB}\rightarrow E}` and
:math:`v_{\mathrm{SSB}\rightarrow E}` are evaluated from planetary
ephemerides in GCRF axes, then

.. math::

   r_\mathrm{ICRF}
   =
   r_\mathrm{GCRF}
   +
   r_{\mathrm{SSB}\rightarrow E},

.. math::

   v_\mathrm{ICRF}
   =
   v_\mathrm{GCRF}
   +
   v_{\mathrm{SSB}\rightarrow E}.

The inverse subtracts the same Earth barycentric state.

Implemented by
``cpp/lupnt/conversions/frame_conversions.cc :: GcrfToIcrf`` /
``IcrfToGcrf`` (a pure translation by the SSB-to-Earth ephemeris state).

GCRF and ITRF
-------------------------------------------------------------------

The transformation from the geocentric celestial frame (GCRF, an Earth-centered
realization of the inertial ICRF axes) to the Earth-fixed ITRF follows the
IAU 2006/2000A CIO-based (Celestial Intermediate Origin) convention.  It is the
product of three rotations, each evaluated at the time scale it physically
depends on:

.. math::

   R_{\mathrm{ITRF},\mathrm{GCRF}}
   =
   \underbrace{R_\mathrm{polar}(t_\mathrm{TT},\,\text{EOP})}_{\text{polar motion}}\;
   \underbrace{R_\mathrm{ERA}(t_\mathrm{UT1})}_{\text{Earth rotation}}\;
   \underbrace{R_\mathrm{PN}(t_\mathrm{TT})}_{\text{precession--nutation}} .

Read right to left, a GCRF vector is first carried to the Celestial Intermediate
Reference System (bias--precession--nutation), then rotated about the Celestial
Intermediate Pole (CIP) through the sidereal Earth Rotation Angle, and finally
tilted into the terrestrial pole by polar motion.

**Precession--nutation** :math:`R_\mathrm{PN}`.  LuPNT reads the CIP unit-vector
coordinates :math:`X,Y` and the CIO locator :math:`s` from tabulated
IAU 2006/2000A (SOFA) series evaluated at TT, and builds the celestial-to-
intermediate matrix directly from them:

.. math::

   R_\mathrm{PN}
   =
   R_z(-s)
   \begin{bmatrix}
     1-aX^2 & -aXY   & -X \\
     -aXY   & 1-aY^2 & -Y \\
     X      & Y      & 1-a\,(X^2+Y^2)
   \end{bmatrix},
   \qquad
   a=\frac{1}{1+\sqrt{1-X^2-Y^2}} .

**Earth rotation** :math:`R_\mathrm{ERA}=R_z(\theta_\mathrm{ERA})` is a rotation
about the CIP by the Earth Rotation Angle, which is *linear* in UT1:

.. math::

   \theta_\mathrm{ERA}(t_\mathrm{UT1})
   =
   2\pi\left(0.7790572732640 + 1.00273781191135448\,D_\mathrm{UT1}\right),

with :math:`D_\mathrm{UT1}` the number of UT1 days since J2000.  UT1 is recovered
from the fitted :math:`\mathrm{UT1}-\mathrm{UTC}` when a SPICE/EOP model covers
``t_tdb``, otherwise from the ``TDB``\ :math:`\to`\ ``UT1`` conversion.  Because
:math:`\theta_\mathrm{ERA}` advances at the Earth rotation rate
:math:`\omega_\oplus\approx7.292115\times10^{-5}\ \mathrm{rad/s}`, this is the
only factor whose time derivative matters for the velocity transform.

**Polar motion** :math:`R_\mathrm{polar}` aligns the CIP with the ITRF pole using
the pole offsets :math:`x_p,y_p` (from the fitted EOP, or the bundled EOP table
sampled at UTC) and the TIO locator
:math:`s' = -47\ \mu''\times(t_\mathrm{TT}/36525\ \text{d})`:

.. math::

   R_\mathrm{polar}
   =
   R_x(-y_p)\,R_y(-x_p)\,R_z(s') .

The position and velocity mapping is then

.. math::

   r_\mathrm{ITRF}
   =
   R_{\mathrm{ITRF},\mathrm{GCRF}} r_\mathrm{GCRF},

.. math::

   v_\mathrm{ITRF}
   =
   R_{\mathrm{ITRF},\mathrm{GCRF}} v_\mathrm{GCRF}
   +
   \dot{R}_{\mathrm{ITRF},\mathrm{GCRF}} r_\mathrm{GCRF},

where the rotation-rate term keeps only the dominant Earth-rotation derivative,

.. math::

   \dot{R}_{\mathrm{ITRF},\mathrm{GCRF}}
   =
   R_\mathrm{polar}\,\dot{R}_\mathrm{ERA}\,R_\mathrm{PN} ,

so :math:`\dot{R}\,r_\mathrm{GCRF}` reproduces the :math:`\omega_\oplus\times r`
transport velocity (the slow precession--nutation and polar-motion rates are
neglected).

The inverse is

.. math::

   r_\mathrm{GCRF}
   =
   R_{\mathrm{ITRF},\mathrm{GCRF}}^\mathsf{T} r_\mathrm{ITRF},

.. math::

   v_\mathrm{GCRF}
   =
   R_{\mathrm{ITRF},\mathrm{GCRF}}^\mathsf{T}
   \left(
     v_\mathrm{ITRF}
     -
     \dot{R}_{\mathrm{ITRF},\mathrm{GCRF}} r_\mathrm{GCRF}
   \right).

Implemented by
``cpp/lupnt/conversions/frame_conversions.cc :: GcrfToItrf`` (and the inverse
``ItrfToGcrf``); the rotation product uses ``RotPolarMotion`` /
``RotSideralMotion`` / ``RotPrecessionNutation``:

.. code-block:: cpp

   Mat3 R_GcrfToItrf     = R_po * R_s     * R_pn;   // R_polar R_ERA R_PN
   Mat3 R_GcrfToItrf_dot = R_po * R_s_dot * R_pn;
   Vec3 r_itrf = R_GcrfToItrf * r_gcrf;
   Vec3 v_itrf = R_GcrfToItrf * v_gcrf + R_GcrfToItrf_dot * r_gcrf;

.. note::

   If ``InitFrameConversionFromSpice`` has fitted a SPICE-derived EOP model
   covering ``t_tdb``, ``RotPolarMotion`` / ``RotSideralMotion`` transparently
   use the fitted EOP so this matches ``spice::ConvertFrameSpice`` closely.

GCRF and EME
-------------------------------------------------------------------

GCRF and EME are related by the Earth frame-bias rotation

.. math::

   R_{\mathrm{EME},\mathrm{GCRF}}
   =
   R_x(-\eta_0) R_y(\xi_0) R_z(\Delta \alpha_0).

The same static rotation is applied to both position and velocity.

Implemented by the static bias matrix
``cpp/lupnt/conversions/frame_conversions.cc :: RotGcrfToEme``, applied by
``GcrfToEme`` / ``EmeToGcrf``:

.. code-block:: cpp

   Mat3d RotGcrfToEme() {
     return RotX(-FRAME_BIAS_ETA0) * RotY(FRAME_BIAS_XI0) * RotZ(FRAME_BIAS_DALPHA0);
   }

GCRF and MOON_CI
-------------------------------------------------------------------

MOON_CI is Moon-centered and axis-aligned with the inertial GCRF axes.  Let
:math:`r_{E\rightarrow M}` and :math:`v_{E\rightarrow M}` be the Moon state
relative to Earth.  Then

.. math::

   r_\mathrm{MOON\_CI}
   =
   r_\mathrm{GCRF}
   -
   r_{E\rightarrow M},

.. math::

   v_\mathrm{MOON\_CI}
   =
   v_\mathrm{GCRF}
   -
   v_{E\rightarrow M}.

The inverse adds the same Earth-to-Moon state.

Implemented by
``cpp/lupnt/conversions/frame_conversions.cc :: GcrfToMoonCi`` /
``MoonCiToGcrf`` (a pure translation by the Earth-to-Moon ephemeris state):

.. code-block:: cpp

   Vec6 GcrfToMoonCi(Real t_tdb, const Vec6& rv_gcrf) {
     Vec6 rv_earth2moon = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
     return rv_gcrf - rv_earth2moon;
   }

MOON_CI and MOON_PA
-------------------------------------------------------------------

MOON_PA is the Moon principal-axis fixed frame.  The native rotation is
represented by lunar mantle libration angles
:math:`(\phi,\theta,\psi)`:

.. math::

   R_{\mathrm{PA},\mathrm{CI}}
   =
   R_z(\psi) R_x(\theta) R_z(\phi).

The velocity transform includes the exact derivative of this rotation:

.. math::

   v_\mathrm{PA}
   =
   R_{\mathrm{PA},\mathrm{CI}} v_\mathrm{CI}
   +
   \dot{R}_{\mathrm{PA},\mathrm{CI}} r_\mathrm{CI}.

The inverse subtracts the rotation-rate contribution before applying the
transpose rotation.

Implemented by
``cpp/lupnt/conversions/frame_conversions.cc :: RotMoonCiToPa`` (the ZXZ angle
reconstruction and its analytic :math:`\dot R`), applied by ``MoonCiToPa`` /
``MoonPaToCi``:

.. code-block:: cpp

   auto [phi, theta, psi, phi_dot, theta_dot, psi_dot] = Unpack(GetLunarMantleData(t_tdb, true));
   Mat3 R_mi2pa = RotZ(psi) * RotX(theta) * RotZ(phi);
   *R_mi2pa_dot = RotZdot(psi, psi_dot) * RotX(theta)   * RotZ(phi)
                + RotZ(psi)   * RotXdot(theta, theta_dot) * RotZ(phi)
                + RotZ(psi)   * RotX(theta)   * RotZdot(phi, phi_dot);

.. note::

   If ``InitFrameConversionFromSpice`` has fitted a lunar-orientation model
   covering ``t_tdb``, the SPICE-fitted MOON_CI->MOON_PA rotation (and its
   exact derivative) is used in place of the DE-Chebyshev libration angles
   (``TryGetFittedMoonCiToPa``).

MOON_PA and MOON_ME
-------------------------------------------------------------------

MOON_ME is related to MOON_PA through a static bias rotation:

.. math::

   R_{\mathrm{ME},\mathrm{PA}}
   =
   R_x(-0.2785'')
   R_y(-78.6944'')
   R_z(-67.8526'').

The same static rotation is applied to both position and velocity.

Implemented by the static matrix
``cpp/lupnt/conversions/frame_conversions.cc :: RotMoonPaToMe``, applied by
``MoonPaToMe`` / ``MoonMeToPa``:

.. code-block:: cpp

   Mat3d RotMoonPaToMe() {
     return RotX(-0.2785 * RAD_ARCSEC) * RotY(-78.6944 * RAD_ARCSEC) * RotZ(-67.8526 * RAD_ARCSEC);
   }

MOON_OP
-------------------------------------------------------------------

MOON_OP is the lunar orbit-plane frame.  It is built from the Moon-to-Earth
state and the Moon mean-Earth pole.  Let

.. math::

   \hat{z}_\mathrm{OP}
   =
   \frac{r_{M\rightarrow E} \times v_{M\rightarrow E}}
        {\lVert r_{M\rightarrow E} \times v_{M\rightarrow E} \rVert},

and let :math:`\hat{p}` be the Moon ME z-axis expressed in MOON_CI.  Then

.. math::

   \hat{x}_\mathrm{OP}
   =
   \frac{\hat{p} \times \hat{z}_\mathrm{OP}}
        {\lVert \hat{p} \times \hat{z}_\mathrm{OP} \rVert},
   \qquad
   \hat{y}_\mathrm{OP}
   =
   \hat{z}_\mathrm{OP} \times \hat{x}_\mathrm{OP}.

The rotation from MOON_OP to MOON_CI is

.. math::

   R_{\mathrm{CI},\mathrm{OP}}
   =
   \begin{bmatrix}
     \hat{x}_\mathrm{OP} &
     \hat{y}_\mathrm{OP} &
     \hat{z}_\mathrm{OP}
   \end{bmatrix}.

The current implementation applies this rotation to both position and velocity
without an explicit :math:`\dot{R}` term for MOON_OP.

Implemented by
``cpp/lupnt/conversions/frame_conversions.cc :: RotMoonOpToCi``, applied by
``MoonOpToCi`` / ``MoonCiToOp``:

.. code-block:: cpp

   Vec6 rv_m2e = GetBodyPosVel(t_tdb, BodyId::MOON, BodyId::EARTH, Frame::GCRF);
   Vec3 z_op   = rv_m2e.head(3).cross(rv_m2e.tail(3)).normalized();
   Vec3 i_pole = ConvertFrame(t_tdb, Vec3::UnitZ(), Frame::MOON_ME, Frame::MOON_CI);
   Vec3 x_op   = i_pole.cross(z_op).normalized();
   Vec3 y_op   = z_op.cross(x_op).normalized();
   R_op2ci << x_op, y_op, z_op;

Solar-System Planet Frames (IAU Orientation)
-------------------------------------------------------------------

Each solar-system ``<PLANET>_CI`` frame is ICRF-aligned and centered on the
planet, and each ``<PLANET>_FIXED`` frame co-rotates with it under the IAU
linear orientation model.  With pole right ascension/declination
:math:`(\alpha, \delta)` and prime-meridian angle :math:`W`, the ICRF-to-fixed
rotation is

.. math::

   R_{\mathrm{FIXED},\mathrm{CI}}
   =
   R_z(W)\,R_x(\tfrac{\pi}{2}-\delta)\,R_z(\tfrac{\pi}{2}+\alpha),

with :math:`(\alpha, \delta)` linear in Julian centuries past J2000 (TT) and
:math:`W` linear in days past J2000 (TT).  Only the prime-meridian spin
contributes appreciably to :math:`\dot{R}`.

Implemented by
``cpp/lupnt/conversions/frame_conversions.cc :: RotBodyCiToFixed`` (with the
per-body coefficients in ``IauOrientationTable``), applied by
``BodyCiToFixed`` / ``BodyFixedToCi``; the frame predicates
``IsPlanetFrame`` / ``IsPlanetFixedFrame`` / ``IsPlanetCiFrame`` classify these
frames, and
``cpp/lupnt/conversions/frame_converter.cc :: ConvertPlanetFrame`` routes them
through the ICRF hub (with a same-body fast path):

.. code-block:: cpp

   // RotBodyCiToFixed
   Real ra  = (o.ra0  + o.ra0_T  * T) * RAD;
   Real dec = (o.dec0 + o.dec0_T * T) * RAD;
   Real W   = (o.w0   + o.w_dot  * d) * RAD;
   Mat3 R = RotZ(W) * RotX(M_PI / 2 - dec) * RotZ(M_PI / 2 + ra);
   if (R_dot) *R_dot = RotZdot(W, o.w_dot * RAD / SECS_DAY) * RotX(M_PI / 2 - dec) * RotZ(M_PI / 2 + ra);

SPICE-Fitted Orientation Option
-------------------------------------------------------------------

``InitFrameConversionFromSpice`` can fit the native Earth and lunar
orientation paths to SPICE over a requested TDB time window.  It samples SPICE
rotations at Chebyshev-Gauss nodes and fits ZXZ Euler angles
:math:`(\phi,\theta,\psi)` over piecewise Chebyshev segments:

.. math::

   q_k(t)
   \approx
   \sum_{n=0}^{N-1} c_{kn} T_n(\tau),
   \qquad
   q_k \in \{\phi,\theta,\psi\}.

At evaluation time the fitted angles reconstruct

.. math::

   R(t) = R_z(\psi(t))R_x(\theta(t))R_z(\phi(t)),

and :math:`\dot{R}(t)` is obtained by differentiating the Chebyshev
polynomials and applying the product rule.  The reconstructed rotation remains
orthonormal to floating-point precision because it is a product of elementary
rotation matrices.

If the requested epoch is outside the fitted time window, LuPNT falls back to
the native analytic Earth or lunar orientation model.

Implemented by
``cpp/lupnt/conversions/frame_conversions.cc :: InitFrameConversionFromSpice``
(the Euler-angle fit is ``FitEulerAngleModel``; the fitted lunar rotation is
consumed by ``TryGetFittedMoonCiToPa``, and the fitted EOP by the Earth
rotation helpers).

Coordinate Scale and Units
-------------------------------------------------------------------

Frame conversion functions themselves do not change numerical units.  If the
input state is in kilometers, the output is in kilometers.  If it is in meters,
the output is in meters.

Ephemeris helpers such as ``GetBodyPosVel`` can additionally scale positions
and gravitational parameters between supported coordinate scales.  For a
constant coordinate scale factor :math:`F`, LuPNT uses

.. math::

   r_\mathrm{to} = F r_\mathrm{from},
   \qquad
   \mu_\mathrm{to} = F \mu_\mathrm{from},
   \qquad
   v_\mathrm{to} = v_\mathrm{from}.

Frame conversion callers are responsible for keeping coordinate scale and unit
choices consistent with the states being transformed.
