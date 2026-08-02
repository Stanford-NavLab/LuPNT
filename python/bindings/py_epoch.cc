// lupnt
#include <lupnt/conversions/epoch.h>
#include <lupnt/conversions/time_conversions.h>
#include <lupnt/core/constants.h>

// pybind11
#include <pybind11/operators.h>

#include <string>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitEpoch(py::module& m) {
  py::class_<Epoch>(m, "Epoch", R"doc(
A time-system-tagged instant stored as exact integer seconds from J2000 plus a
fractional second, so differences and time-scale conversions keep sub-nanosecond
precision at any date.

A bare float epoch (seconds from J2000) has a ULP of |t|*2**-52 -- about 0.25 us
at present-day dates, i.e. ~75 m at the speed of light. Any small quantity pulled
back out of such a value is snapped to that grid. `Epoch` avoids this by never
storing the large and small parts in the same float.

    e  = pnt.Epoch.from_gregorian("2030-01-01T00:00:00", pnt.Time.TAI)
    e2 = e + 1e-12                 # 1 ps later -- representable
    (e2 - e)                       # -> 1e-12 exactly; a float epoch gives 0.0
    e.to(pnt.Time.TDB)             # scale conversion, sub-ns

`to_seconds()` collapses back to a bare float for the existing API surface; that
step is lossy by construction, so prefer `-` and `to()` when precision matters.
)doc")
      .def(py::init<int64_t, Real, Time>(), py::arg("sec"), py::arg("frac"), py::arg("scale"),
           "Construct from integer seconds + fractional second (need not be normalised).")
      .def_static("from_seconds", &Epoch::FromSeconds, py::arg("t"), py::arg("scale"),
                  "Build from seconds-from-J2000. Lossy: cannot recover bits `t` never had.")
      .def_static(
          "from_gregorian",
          [](int y, int mo, int d, int h, int mi, Real s, Time scale) {
            return Epoch::FromGregorian(y, mo, d, h, mi, s, scale);
          },
          py::arg("year"), py::arg("month"), py::arg("day"), py::arg("hour") = 0,
          py::arg("minute") = 0, py::arg("second") = Real(0.0), py::arg("scale") = Time::TAI,
          "Build from calendar fields; the seconds field may be fractional.")
      .def_static(
          "from_gregorian",
          [](const std::string& date, Time scale) { return Epoch::FromGregorian(date, scale); },
          py::arg("date"), py::arg("scale") = Time::TAI,
          "Build from an ISO-like string, e.g. '2030-01-01T00:00:00'.")

      .def_property_readonly("sec", &Epoch::seconds, "Exact integer seconds from J2000.")
      .def_property_readonly("frac", &Epoch::fraction, "Fractional second in [0, 1).")
      .def_property_readonly("scale", &Epoch::scale, "Time scale this epoch is expressed in.")

      .def("to_seconds", &Epoch::ToSeconds,
           "Collapse to seconds from J2000. Lossy -- reintroduces the ~0.25 us epoch ULP.")
      .def("to_mjd", &Epoch::ToMjd, "Modified Julian Date in this epoch's scale (lossy).")
      .def("to", &Epoch::To, py::arg("scale"),
           "Convert to another time system, retaining sub-nanosecond precision.")
      .def("to_gregorian_string", &Epoch::ToGregorianString, py::arg("precision") = 6,
           "Calendar string in this epoch's own scale.")

      // NOTE: register `Epoch - Epoch` BEFORE `Epoch - Real`. pybind11 tries
      // overloads in registration order, and LuPNT's `Real` caster raises
      // TypeError (rather than returning false) when handed a non-numeric
      // object -- so if the Real overload were tried first it would poison the
      // dispatch with a live exception instead of falling through.
      .def(
          "__sub__", [](const Epoch& a, const Epoch& b) { return a - b; }, py::is_operator(),
          "Elapsed time [s] between two epochs in the same scale.")
      .def(
          "__sub__", [](const Epoch& a, Real dt) { return a - dt; }, py::is_operator(),
          "Retreat by a duration [s].")
      .def(
          "__add__", [](const Epoch& a, Real dt) { return a + dt; }, py::is_operator(),
          "Advance by a duration [s].")
      .def(
          "__radd__", [](const Epoch& a, Real dt) { return a + dt; }, py::is_operator(),
          "Advance by a duration [s].")
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self > py::self)
      .def(py::self <= py::self)
      .def(py::self >= py::self)

      .def("__repr__", [](const Epoch& e) {
        return "<Epoch " + e.ToGregorianString(9) + " " + time_to_string.at(e.scale()) + ">";
      });

  m.def("time_scale_offset", &TimeScaleOffset, py::arg("epoch"), py::arg("to"),
        "Offset (to - from) [s] at this instant, as a small full-precision quantity.");

  py::class_<EpochSeries>(m, "EpochSeries", R"doc(
A sequence of epochs sharing one time scale, stored columnwise as exact integer
seconds plus a vector of fractions.

Use `since(ref)` to feed vectorized APIs without precision loss: it removes the
large common epoch once, so the callee receives small offsets.
)doc")
      .def_static("linspace", &EpochSeries::Linspace, py::arg("start"), py::arg("step"),
                  py::arg("n"),
                  "`n` epochs spaced `step` seconds apart; accumulated in integer seconds so a "
                  "long grid does not drift.")
      .def_static("from_seconds", &EpochSeries::FromSeconds, py::arg("t"), py::arg("scale"),
                  "Wrap a bare seconds-from-J2000 array (lossy, for interop).")
      .def("__len__", &EpochSeries::size)
      .def("__getitem__", &EpochSeries::operator[], py::arg("i"))
      .def_property_readonly("scale", &EpochSeries::scale)
      .def("to_seconds", &EpochSeries::ToSeconds, "Collapse to bare seconds from J2000 (lossy).")
      .def("since", &EpochSeries::Since, py::arg("ref"),
           "Elapsed seconds of every sample relative to `ref`, exactly.")
      .def("to", &EpochSeries::To, py::arg("scale"), "Convert the whole series to another scale.");
}
