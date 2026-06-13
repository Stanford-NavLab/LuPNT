#include "lupnt/environment/forces.h"

#include <omp.h>

#include "lupnt/conversions/coordinate_conversions.h"
#include "lupnt/conversions/time_conversions.h"
#include "lupnt/environment/solar_system.h"

namespace lupnt {
  namespace {
    Real ClampUnit(Real x) {
      if (x.val() > 1.0) return 1.0;
      if (x.val() < -1.0) return -1.0;
      return x;
    }

    Real Clamp01(Real x) {
      if (x.val() > 1.0) return 1.0;
      if (x.val() < 0.0) return 0.0;
      return x;
    }
  }  // namespace

  /// @brief Computes the acceleration due to the harmonic gravity field of the
  /// central body
  /// @param r Satellite position in the body fixed system [m]
  /// @param GM Gravitational coefficient [m^3/s^2]
  /// @param R_ref Reference radius [m]
  /// @param CS Spherical harmonics coefficients (un-normalized)
  /// @param n_max Maximum degree
  /// @param m_max Maximum order (m_max<=n_max; m_max=0 for zonals, only)
  /// @return Acceleration [m/s^2]
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  template <typename T> Vector3<T> AccelarationGravityField(const Vector3<T>& r, T GM, T R_ref,
                                                            const MatrixX<T>& CS, int n_max,
                                                            int m_max) {
    MatrixX<T> V(n_max + 2, n_max + 2);  // Harmonic functions
    MatrixX<T> W(n_max + 2, n_max + 2);  // work array (0..n_max+1,0..n_max+1)

    T r_sqr = r.squaredNorm();      // [m^2]
    T rho = R_ref * R_ref / r_sqr;  // [-]

    // Normalized coordinates
    T x0 = R_ref * r(0) / r_sqr;  // [-]
    T y0 = R_ref * r(1) / r_sqr;  // [-]
    T z0 = R_ref * r(2) / r_sqr;  // [-]

    // Harmonic functions up to degree and order n_max+1
    //   V_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
    //   W_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)

    // Zonal terms V(n,0); set W(n,0)=0.0
    V(0, 0) = R_ref / sqrt(r_sqr);  // [-]
    V(1, 0) = z0 * V(0, 0);         // [-]
    W(0, 0) = 0.0;                  // [-]
    W(1, 0) = 0.0;                  // [-]

    for (int n = 2; n <= n_max + 1; n++) {
      V(n, 0) = ((2 * n - 1) * z0 * V(n - 1, 0) - (n - 1) * rho * V(n - 2, 0)) / n;
      W(n, 0) = 0.0;
    };

    // Tesseral and sectorial terms
    for (int m = 1; m <= m_max + 1; m++) {
      // V(m,m) .. V(n_max+1,m)
      V(m, m) = (2 * m - 1) * (x0 * V(m - 1, m - 1) - y0 * W(m - 1, m - 1));
      W(m, m) = (2 * m - 1) * (x0 * W(m - 1, m - 1) + y0 * V(m - 1, m - 1));
      if (m <= n_max) {
        V(m + 1, m) = (2 * m + 1) * z0 * V(m, m);
        W(m + 1, m) = (2 * m + 1) * z0 * W(m, m);
      };

      for (int n = m + 2; n <= n_max + 1; n++) {
        V(n, m) = ((2 * n - 1) * z0 * V(n - 1, m) - (n + m - 1) * rho * V(n - 2, m)) / (n - m);
        W(n, m) = ((2 * n - 1) * z0 * W(n - 1, m) - (n + m - 1) * rho * W(n - 2, m)) / (n - m);
      };
    };

    T ax = 0, ay = 0, az = 0;
    // #pragma omp parallel for reduction(+ : ax, ay, az)
    for (int m = 0; m <= m_max; m++)
      for (int n = m; n <= n_max; n++)
        if (m == 0) {
          T C = CS(n, 0);  // C_n,0
          ax -= C * V(n + 1, 1);
          ay -= C * W(n + 1, 1);
          az -= (n + 1) * C * V(n + 1, 0);
        } else {
          T C = CS(n, m);      // C_n,m
          T S = CS(m - 1, n);  // S_n,m
          T Fac = 0.5 * (n - m + 1) * (n - m + 2);
          ax += +0.5 * (-C * V(n + 1, m + 1) - S * W(n + 1, m + 1))
                + Fac * (+C * V(n + 1, m - 1) + S * W(n + 1, m - 1));
          ay += +0.5 * (-C * W(n + 1, m + 1) + S * V(n + 1, m + 1))
                + Fac * (-C * W(n + 1, m - 1) + S * V(n + 1, m - 1));
          az += (n - m + 1) * (-C * V(n + 1, m) - S * W(n + 1, m));
        };

    Vector3<T> a = (GM / (R_ref * R_ref)) * Vector3<T>(ax, ay, az);
    return a;
  }

  template Vec3d AccelarationGravityField(const Vec3d& r, double GM, double R_ref, const MatXd& CS,
                                          int n_max, int m_max);
  template Vec3 AccelarationGravityField(const Vec3& r, Real GM, Real R_ref, const MatX& CS,
                                         int n_max, int m_max);

  /// @brief Computes the acceleration due to a point mass
  /// @param r Satellite position [m]
  /// @param s Point mass position [m]
  /// @param GM Gravitational coefficient of point mass [m^3/s^2]
  /// @return Acceleration [m/s^2]
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Vec3 AccelerationPointMass(const Vec3& r, const Vec3& s, Real GM) {
    Vec3 d = r - s;
    // a = -GM * (d / |d|^3 + s / |s|^3)
    Vec3 a = Vec3::Zero();
    if (s.norm() > EPS) a += s / pow(s.norm(), 3);
    if (d.norm() > EPS) a += d / pow(d.norm(), 3);
    a *= -GM;
    return a;
  }

  /// @brief Computes the acceleration due to solar radiation pressure
  /// @param r Spacecraft position [m]
  /// @param r_sun Sun position [m]
  /// @param area Cross-section [m^2]
  /// @param mass Spacecraft mass [kg]
  /// @param CR Solar radiation pressure coefficient (0: translucent, 1: black
  /// body 2: perfect mirror)
  /// @param P0 Solar radiation pressure at 1 AUa
  /// @param AU Length of one Astronomical Unit
  /// @return Acceleration [m/s^2]
  /// @note
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Vec3 AccelerationSolarRadiation(const Vec3& r, const Vec3& r_sun, Real bcoeff_srp, Real P0,
                                  Real AU) {
    Vec3 d = r - r_sun;

    // Vec3 a = CR * (area / mass) * P0 * (AU * AU) * d / pow(d.norm(), 3);
    Vec3 a = bcoeff_srp * P0 * (AU * AU) * d / pow(d.norm(), 3);

    return a;
  }

  /// @brief Computes the first-order relativistic correction due to a central body
  /// @param r Spacecraft position relative to central body [m]
  /// @param v Spacecraft velocity relative to central body [m/s]
  /// @param GM Gravitational coefficient of central body [m^3/s^2]
  /// @param c_light Speed of light [m/s]
  /// @return Relativistic acceleration correction [m/s^2]
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications, Sec. 3.7.3 Eq. 3.146. Berlin: New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Vec3 AccelerationRelativisticCorrection(const Vec3& r, const Vec3& v, Real GM, Real c_light) {
    Real r_norm = r.norm();
    if (r_norm <= EPS) return Vec3::Zero();

    Real c2 = c_light * c_light;
    Real v2 = v.squaredNorm();
    Real rv = r.dot(v);
    return -GM / (c2 * pow(r_norm, 3)) * ((4.0 * GM / r_norm - v2) * r + 4.0 * rv * v);
  }

  /// @brief Computes the acceleration due to the atmospheric drag
  /// @param mjd_tt Terrestrial Time (Modified Julian Date)
  /// @param rv Satellite position and velocity in the inertial system [km, km/s]
  /// @param R Transformation matrix to true-of-date inertial system
  /// @param area Cross-section [m^2]
  /// @param mass Spacecraft mass [kg]
  /// @param CD Drag coefficient
  /// @return Acceleration [km/s^2]
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Vec3 AccelerationDrag(Real mjd_tt, const Vec6& rv, const Mat3& T, Real bcoeff_drag) {
    const Vec3 omega(0, 0, 7.29212e-5);  // Earth angular velocity [rad/s]
    Vec3 r = rv.head(3);
    Vec3 v = rv.tail(3);

    // Position and velocity in true-of-date system
    Vec3 r_tod = T * r;
    Vec3 v_tod = T * v;

    // Velocity relative to the Earth's atmosphere
    Vec3 v_rel = v_tod - omega.cross(r_tod);
    Real v_abs = v_rel.norm();

    // Atmospheric density due to modified Harris-Priester model
    Real dens = DensityHarrisPriester(mjd_tt, r_tod);

    // Acceleration
    // Vec3 a_tod = -0.5 * CD * (area / mass) * dens * v_abs * v_rel * KM_M;
    Vec3 a_tod = -0.5 * bcoeff_drag * dens * v_abs * v_rel * KM_M;

    // Transformation  to ICRF/EME2000 system
    return T.transpose() * a_tod;
  }

  /// @brief Computes the atmospheric density for the modified Harris-Priester
  /// model
  /// @param mjd_tt Terrestrial Time (Modified Julian Date)
  /// @param r_tod Satellite position in the inertial system [m]
  /// @return Density [kg/m^3]
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Real DensityHarrisPriester(Real mjd_tt, const Vec3& r_tod) {
    const double upper_limit = 1000.0e3;  // Upper height limit [m]
    const double lower_limit = 100.0e3;   // Lower height limit [m]
    const double ra_lag = 0.523599;       // Right ascension lag [rad]
    const int n_prm = 3;                  // Harris-Priester parameter
                                          // 2(6) low(high) inclination

    // Harris-Priester atmospheric density model parameters
    // Height [m], minimum density, maximum density [gm/km^3]
    const int N_Coef = 50;
    const Vecd<N_Coef> h = {100e3, 120e3, 130e3, 140e3, 150e3, 160e3, 170e3, 180e3, 190e3, 200e3,
                            210e3, 220e3, 230e3, 240e3, 250e3, 260e3, 270e3, 280e3, 290e3, 300e3,
                            320e3, 340e3, 360e3, 380e3, 400e3, 420e3, 440e3, 460e3, 480e3, 500e3,
                            520e3, 540e3, 560e3, 580e3, 600e3, 620e3, 640e3, 660e3, 680e3, 700e3,
                            720e3, 740e3, 760e3, 780e3, 800e3, 840e3, 880e3, 920e3, 960e3, 1000e3};
    const Vecd<N_Coef> c_min
        = {4.974e+05, 2.490e+04, 8.377e+03, 3.899e+03, 2.122e+03, 1.263e+03, 8.008e+02, 5.283e+02,
           3.617e+02, 2.557e+02, 1.839e+02, 1.341e+02, 9.949e+01, 7.488e+01, 5.709e+01, 4.403e+01,
           3.430e+01, 2.697e+01, 2.139e+01, 1.708e+01, 1.099e+01, 7.214e+00, 4.824e+00, 3.274e+00,
           2.249e+00, 1.558e+00, 1.091e+00, 7.701e-01, 5.474e-01, 3.916e-01, 2.819e-01, 2.042e-01,
           1.488e-01, 1.092e-01, 8.070e-02, 6.012e-02, 4.519e-02, 3.430e-02, 2.632e-02, 2.043e-02,
           1.607e-02, 1.281e-02, 1.036e-02, 8.496e-03, 7.069e-03, 4.680e-03, 3.200e-03, 2.210e-03,
           1.560e-03, 1.150e-03};
    const Vecd<N_Coef> c_max
        = {4.974e+05, 2.490e+04, 8.710e+03, 4.059e+03, 2.215e+03, 1.344e+03, 8.758e+02, 6.010e+02,
           4.297e+02, 3.162e+02, 2.396e+02, 1.853e+02, 1.455e+02, 1.157e+02, 9.308e+01, 7.555e+01,
           6.182e+01, 5.095e+01, 4.226e+01, 3.526e+01, 2.511e+01, 1.819e+01, 1.337e+01, 9.955e+00,
           7.492e+00, 5.684e+00, 4.355e+00, 3.362e+00, 2.612e+00, 2.042e+00, 1.605e+00, 1.267e+00,
           1.005e+00, 7.997e-01, 6.390e-01, 5.123e-01, 4.121e-01, 3.325e-01, 2.691e-01, 2.185e-01,
           1.779e-01, 1.452e-01, 1.190e-01, 9.776e-02, 8.059e-02, 5.741e-02, 4.210e-02, 3.130e-02,
           2.360e-02, 1.810e-02};

    // Satellite height (Earth flattening correction)
    Vec3 lla = CartToLatLonAlt(r_tod, R_EARTH, WGS84_F);
    Real height = lla(2);  // [m]

    // Exit with zero density outside height model limits
    if (height >= upper_limit || height <= lower_limit) {
      return 0.0;
    }

    // Sun right ascension, declination
    Vec3 r_sun = SunPositionLowPrecision(mjd_tt);
    Real ra_sun = atan2(r_sun(1), r_sun(0));
    Real dec_sun = atan2(r_sun(2), sqrt(pow(r_sun(0), 2) + pow(r_sun(1), 2)));

    // Unit vector u towards the apex of the diurnal bulge
    // in inertial geocentric coordinates
    Real c_dec = cos(dec_sun);
    Vec3 u;  // Apex of diurnal bulge
    u(0) = c_dec * cos(ra_sun + ra_lag);
    u(1) = c_dec * sin(ra_sun + ra_lag);
    u(2) = sin(dec_sun);

    // Cosine of half angle between satellite position and
    // apex of diurnal bulge (Harris-Priester modification )
    Real c_psi2 = 0.5 + 0.5 * r_tod.dot(u) / r_tod.norm();

    // Height index search and exponential density interpolation
    int ih = 0;  // section index reset
    // loop over N_Coef height regimes
    for (int i = 0; i < N_Coef - 1; i++) {
      if (height >= h(i) && height < h(i + 1)) {
        ih = i;  // ih identifies height section
        break;
      }
    }

    Real h_min = (h(ih) - h(ih + 1)) / log(c_min(ih + 1) / c_min(ih));
    Real h_max = (h(ih) - h(ih + 1)) / log(c_max(ih + 1) / c_max(ih));

    Real d_min = c_min(ih) * exp((h(ih) - height) / h_min);
    Real d_max = c_max(ih) * exp((h(ih) - height) / h_max);

    // Density computation
    Real density = d_min + (d_max - d_min) * pow(c_psi2, n_prm);
    return density * 1.0e-9;  // [kg/m^3]
  }

  /// @brief Computes the solar-radiation shadow function of a spacecraft
  /// @param r Spacecraft position relative to occulting body [m]
  /// @param r_sun Sun position relative to occulting body [m]
  /// @param R_body Radius of occulting body [m]
  /// @param R_sun Radius of Sun [m]
  /// @return Shadow function [0,1], 0 if in umbra, 1 if fully illuminated
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications, Sec. 3.4.2. Berlin: New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Real ShadowFunction(const Vec3& r, const Vec3& r_sun, Real R_body, Real R_sun) {
    Real r_norm = r.norm();
    Vec3 rho_sun = r_sun - r;  // spacecraft -> Sun
    Real rho_sun_norm = rho_sun.norm();
    if (r_norm <= EPS || rho_sun_norm <= EPS) return 1.0;

    Real a = asin(ClampUnit(R_sun / rho_sun_norm));  // apparent solar radius
    Real b = asin(ClampUnit(R_body / r_norm));       // apparent occulting-body radius
    Real c = acos(ClampUnit(-r.dot(rho_sun) / (r_norm * rho_sun_norm)));

    if (c >= a + b) return 1.0;

    if (c <= abs(b - a)) {
      if (b >= a) return 0.0;
      return Clamp01(1.0 - b * b / (a * a));
    }

    Real x = (c * c + a * a - b * b) / (2.0 * c);
    Real y2 = a * a - x * x;
    Real y = sqrt(Max(y2, 0.0));

    Real occulted_area
        = a * a * acos(ClampUnit(x / a)) + b * b * acos(ClampUnit((c - x) / b)) - c * y;
    return Clamp01(1.0 - occulted_area / (PI * a * a));
  }

  /// @brief Computes the fractional illumination of a spacecraft
  /// @param r Spacecraft position relative to occulting body [m]
  /// @param r_sun Sun position relative to occulting body [m]
  /// @return Illumination factor [0,1], 0 if in shadow, 1 if fully illuminated
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications, Sec. 3.4.2. Berlin: New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Real Illumination(const Vec3& r, const Vec3& r_sun, Real R_body) {
    return ShadowFunction(r, r_sun, R_body, R_SUN);
  }

  /// @brief Computes the acceleration acting on a satellite orbiting the Earth
  /// @param mjd_tt Terrestrial Time (Modified Julian Date)
  /// @param rv Satellite position and velocity in the inertial system [m]
  /// @param area Cross-section [m^2]
  /// @param mass Spacecraft mass [kg]
  /// @param CR Solar radiation pressure coefficient [N/m^2]
  /// @param CD Drag coefficient [-]
  /// @return Acceleration [m/s^2]
  /// @return Illumination factor [0,1], 0 if in shadow, 1 if fully illuminated
  /// @ref
  /// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
  /// applications. Berlin : New York: Springer, 2000.
  /// doi: 10.1007/978-3-642-58351-3.
  Vec3 AccelerationEarthSpacecraft(Real mjd_tt, const Vec6& rv, Real bcoeff_srp, Real bcoeff_drag,
                                   GravityField<Real> grav) {
    Vec3 r = rv.head(3);

    // Acceleration due to harmonic gravity field
    Real mjd_ut1 = mjd_tt;
    Mat3 T = NutationMatrix(mjd_tt) * PrecessionMatrix(MJD_J2000_TT, mjd_tt);
    Mat3 E = GreenwichHourAngleMatrix(mjd_ut1) * T;
    Vec3 r_bf = E * r;
    Vec3 a_bf = AccelarationGravityField(r_bf, grav.GM, grav.R, grav.CS, grav.n, grav.m);
    Vec3 a = T.transpose() * a_bf;

    // Luni-solar perturbations
    Vec3 r_sun = SunPositionLowPrecision(mjd_tt);
    Vec3 r_Moon = MoonPositionLowPrecision(mjd_tt);
    a += AccelerationPointMass(r, r_sun, GM_SUN);
    a += AccelerationPointMass(r, r_Moon, GM_MOON);

    // Relativistic correction from the central Earth mass
    a += AccelerationRelativisticCorrection(r, rv.tail(3), grav.GM);

    // Solar radiation pressure
    a += Illumination(r, r_sun, R_EARTH)
         * AccelerationSolarRadiation(r, r_sun, bcoeff_srp, P_SUN, AU);

    // Atmospheric drag
    a += AccelerationDrag(mjd_tt, rv, T, bcoeff_drag);

    return a;
  }

}  // namespace lupnt
