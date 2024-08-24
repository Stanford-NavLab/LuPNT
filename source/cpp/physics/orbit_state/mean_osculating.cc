#include "lupnt/physics/orbit_state/mean_osculating.h"

#include "lupnt/physics/orbit_state/conversions.h"

namespace lupnt {

  // Mean and Osculating
  Vec6 Mean2Osculating(const Vec6 &coe_m, Real GM, Real J2) {
    Vec6 coe_o;

    if (J2 > 0) {
      Vec6 meanEquioe = Classical2Equinoctial(coe_m, GM);
      Vec6 oscEquioe;  // = MeanOscClosedEqui(meanEquioe, J2);
      coe_o = Equinoctial2Classical(oscEquioe, GM);
    } else {
      coe_o = coe_m;
    }

    return coe_o;
  }

  ClassicalOE Mean2Osculating(const ClassicalOE &coe_m, Real GM, Real J2) {
    return ClassicalOE(Mean2Osculating(coe_m.GetVec(), GM, J2), coe_m.GetFrame());
  }

  Vec6 osc2mean_NRiterator(const Vec6 &osc_equi_elem, double tol) {
    Vec6 mean_equi_elem = osc_equi_elem;
    double R = 1.0;
    int niter = 0;

    while (std::abs(R) > tol) {
      niter++;
      Vec6 osc_loop;
      // std::tie(std::ignore, osc_loop, std::ignore) =
      // transformationmatrix_osc2mean_equinoctial(mean_equi_elem);
      Vec6 delta = osc_equi_elem - osc_loop;
      // R = delta.norm_inf(); // Assuming the library provides an infinity
      // norm function
      mean_equi_elem = mean_equi_elem + delta;

      if (niter > 100) {
        std::cout << "Osc2Mean iterations > 100" << std::endl;
        break;
      }
    }

    return mean_equi_elem;
  }

  Vec6 Osculating2Mean(const Vec6 &coe_o, Real GM, Real J2) {
    Vec6 coe_m;
    double tol = 1e-8;

    if (J2 > 0) {
      Vec6 eqoe_o = Classical2Equinoctial(coe_o, GM);
      Vec6 eqoe_m = osc2mean_NRiterator(eqoe_o, tol);
      coe_m = Equinoctial2Classical(eqoe_m, GM);
    } else {
      coe_m = coe_o;
    }

    return coe_m;
  }

  ClassicalOE Osculating2Mean(const ClassicalOE &coe_o, Real GM, Real J2) {
    return ClassicalOE(Osculating2Mean(coe_o.GetVec(), GM, J2), coe_o.GetFrame());
  }

}  // namespace lupnt
