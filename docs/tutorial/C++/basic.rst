.. _cplusplus_basic_tutorial:

C++ ELFO Propagation and Frame Conversion
=========================================

This example propagates a spacecraft in a lunar elliptical frozen orbit (ELFO).
The initial state is defined with classical orbital elements in the Moon orbital
plane frame, converted to the Moon-centered inertial frame for dynamics, and then
converted to the Moon principal-axes frame for inspection.

.. code-block:: cpp

    #include <iostream>

    #include "lupnt/lupnt.h"

    using namespace lupnt;

    int main() {
      const Real t0_tdb =
          ConvertTime(GregorianToTime(2030, 1, 1, 12, 0, 0), Time::UTC, Time::TDB);

      Vec6 coe_op;
      coe_op << 6541.4e3, 0.60, 65.5 * RAD, 0.0 * RAD, 90.0 * RAD, 0.0 * RAD;

      const Vec6 rv0_op = ClassicalToCart(coe_op, GM_MOON);
      const Vec6 rv0_ci = ConvertFrame(t0_tdb, rv0_op, Frame::MOON_OP, Frame::MOON_CI);

      NBodyDynamics dynamics;
      dynamics.SetFrame(Frame::MOON_CI);
      dynamics.AddBody(Body::Moon(20, 20));
      dynamics.AddBody(Body::Earth());
      dynamics.AddBody(Body::Sun());

      const VecX tspan = Arange(t0_tdb, t0_tdb + 7.0 * SECS_DAY, SECS_HOUR);
      const MatX rv_ci = dynamics.Propagate(rv0_ci, tspan);
      const MatX rv_pa = ConvertFrame(tspan, rv_ci, Frame::MOON_CI, Frame::MOON_PA);

      std::cout << "Initial Moon-CI state [m, m/s]:\n" << rv0_ci.transpose() << "\n\n";
      std::cout << "Final Moon-PA state [m, m/s]:\n" << rv_pa.row(rv_pa.rows() - 1) << "\n";

      return 0;
    }

To compile this as a standalone example, include the LuPNT headers and link
against the ``lupnt`` library target from the CMake build.
