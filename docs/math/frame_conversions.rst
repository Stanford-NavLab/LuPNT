Frame Conversions
===================================================================

Purpose
-------------------------------------------------------------------

This specification defines the mathematical contract for LuPNT frame
conversion functions such as ``ConvertFrame``,
``GetFrameRotationTranslation``, and the Earth/Moon frame-specific conversion
helpers.  The key requirements are explicit epoch scale, frame origin,
rotation, translation, and velocity conventions.

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

GCRF and ITRF
-------------------------------------------------------------------

The native GCRF-to-ITRF rotation is

.. math::

   R_{\mathrm{ITRF},\mathrm{GCRF}}
   =
   R_\mathrm{polar}
   R_\mathrm{ERA}
   R_\mathrm{PN}.

Here :math:`R_\mathrm{PN}` is the IAU precession-nutation rotation,
:math:`R_\mathrm{ERA}` is Earth rotation from UT1, and
:math:`R_\mathrm{polar}` is polar motion from EOP data.

The position and velocity mapping is

.. math::

   r_\mathrm{ITRF}
   =
   R_{\mathrm{ITRF},\mathrm{GCRF}} r_\mathrm{GCRF},

.. math::

   v_\mathrm{ITRF}
   =
   R_{\mathrm{ITRF},\mathrm{GCRF}} v_\mathrm{GCRF}
   +
   \dot{R}_{\mathrm{ITRF},\mathrm{GCRF}} r_\mathrm{GCRF}.

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

GCRF and EME
-------------------------------------------------------------------

GCRF and EME are related by the Earth frame-bias rotation

.. math::

   R_{\mathrm{EME},\mathrm{GCRF}}
   =
   R_x(-\eta_0) R_y(\xi_0) R_z(\Delta \alpha_0).

The same static rotation is applied to both position and velocity.

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
