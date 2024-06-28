#pragma once
#include "lupnt/numerics/vector_macros.h"
#include "orbit_states.h"

namespace lupnt {
// From CartesianOrbitState
// - 2 ClassicalOE
ClassicalOE Cart2Classical(const CartesianOrbitState &rv, Real GM);
Vec6 Cart2Classical(const Vec6 &rv, Real GM);
Vec6 Cart2Classical(Real dt, const Vec3 &r1, const Vec3 &r2, Real GM);

// - 2 CartesianOrbitState (relative)
CartesianOrbitState Inertial2Rtn(const CartesianOrbitState &rv_c,
                                 const CartesianOrbitState &rv_d);
Vec6 Inertial2Rtn(const Vec6 &rv_c, const Vec6 &rv_d);

CartesianOrbitState Rtn2Inertial(const CartesianOrbitState &rv_c,
                                 const CartesianOrbitState &rv_rtn_d);
Vec6 Rtn2Inertial(const Vec6 &rv_c, const Vec6 &rv_rtn_d);

// From ClassicalOE
// - 2 CartesianOrbitState
CartesianOrbitState Classical2Cart(const ClassicalOE &coe, Real GM);
Vec6 Classical2Cart(const Vec6 &coe, Real GM);

// - 2 QuasiNonsingOE
QuasiNonsingOE Classical2QuasiNonsing(const ClassicalOE &coe, Real GM);
Vec6 Classical2QuasiNonsing(const Vec6 &coe, Real GM);

// - 2 EquinoctialOE
EquinoctialOE Classical2Equinoctial(const ClassicalOE &coe, Real GM);
Vec6 Classical2Equinoctial(const Vec6 &coe, Real GM);

// - 2 DelaunayOE
DelaunayOE Classical2Delaunay(const ClassicalOE &coe, Real GM);
Vec6 Classical2Delaunay(const Vec6 &coe, Real GM);

// From QuasiNonsingOE
// - 2 ClassicalOE
ClassicalOE QuasiNonsing2Classical(const QuasiNonsingOE &qnsoe, Real GM);
Vec6 QuasiNonsing2Classical(const Vec6 &qnsoeVec, Real GM);

// From EquinoctialOE
// - 2 ClassicalOE
ClassicalOE Equinoctial2Classical(const EquinoctialOE &eqoe, Real GM);
Vec6 Equinoctial2Classical(const Vec6 &eqoe, Real GM);

// From DelaunayOE
// - 2 ClassicalOE
ClassicalOE Delaunay2Classical(const DelaunayOE &deloe, Real GM);
Vec6 Delaunay2Classical(const Vec6 &deloe, Real GM);

// From ClassicalOE and QuasiNonsingROE (relative)
// - 2 ClassicalOE
ClassicalOE RelQuasiNonsing2Classical(const ClassicalOE &coe,
                                      const QuasiNonsingROE &RelQuasiNonsing);
Vec6 RelQuasiNonsing2Classical(const Vec6 &coe, const Vec6 &RelQuasiNonsing);

// Vector definitions
VEC_DEF_VECTOR_REAL(Classical2Cart, 6);
VEC_DEF_VECTOR_VECTOR(Inertial2Rtn, 6);
VEC_DEF_VECTOR_VECTOR(Rtn2Inertial, 6);
VEC_DEF_VECTOR_REAL(Cart2Classical, 6);
VEC_DEF_VECTOR_REAL(Classical2QuasiNonsing, 6);
VEC_DEF_VECTOR_REAL(Classical2Equinoctial, 6);
VEC_DEF_VECTOR_REAL(Classical2Delaunay, 6);
VEC_DEF_VECTOR_REAL(QuasiNonsing2Classical, 6);
VEC_DEF_VECTOR_REAL(Equinoctial2Classical, 6);
VEC_DEF_VECTOR_REAL(Delaunay2Classical, 6);
VEC_DEF_VECTOR_VECTOR(RelQuasiNonsing2Classical, 6);

}  // namespace lupnt