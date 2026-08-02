// Epoch must build a calendar date exactly.
//
// A calendar instant at a whole second is an exact integer number of seconds
// from J2000, and Epoch stores the integer part as int64 -- so there is no
// reason for any rounding at all. Routing through a Modified Julian Date would
// introduce some: MJD counts days from 1858, so one ULP of an MJD at
// present-day epochs is ~1.1 us, and that rounding lands on the largest
// quantity in the chain. Whole minutes survive it (dyadic fraction-of-day) but
// arbitrary clock times do not.

#include <lupnt/conversions/epoch.h>
#include <lupnt/conversions/time_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Reference: exact seconds from J2000 (2000-01-01 12:00) by integer arithmetic.
  int64_t ExactSeconds(int64_t y, int64_t mo, int64_t d, int64_t h, int64_t mi, int64_t s) {
    int64_t yy = y - (mo <= 2);
    int64_t era = (yy >= 0 ? yy : yy - 399) / 400;
    int64_t yoe = yy - era * 400;
    int64_t doy = (153 * (mo + (mo > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    int64_t doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    int64_t days = era * 146097 + doe - 719468 - 10957;  // days from 2000-01-01
    return days * 86400 + h * 3600 + mi * 60 + s - 43200;
  }
}  // namespace

TEST_CASE("conversions.epoch_calendar_is_exact") {
  struct Case {
    int y, mo, d, h, mi;
    double s;
  };
  // Deliberately includes non-dyadic clock times, which are the ones an
  // MJD-routed conversion gets wrong.
  const Case cases[] = {
      {2020, 1, 1, 12, 0, 0.0},    {2020, 1, 1, 0, 0, 0.0},   {2020, 7, 4, 12, 34, 56.0},
      {2030, 3, 15, 6, 0, 0.0},    {2025, 11, 7, 3, 17, 9.0}, {1999, 12, 31, 23, 59, 59.0},
      {2049, 2, 28, 13, 46, 51.0}, {2000, 1, 1, 12, 0, 0.0},
  };

  for (const Case& c : cases) {
    INFO(c.y << "-" << c.mo << "-" << c.d << "T" << c.h << ":" << c.mi << ":" << c.s);
    Epoch e = Epoch::FromGregorian(c.y, c.mo, c.d, c.h, c.mi, Real(c.s), Time::TDB);
    const int64_t ref = ExactSeconds(c.y, c.mo, c.d, c.h, c.mi, static_cast<int64_t>(c.s));
    // The integer second must be reproduced BIT-EXACTLY, not to a tolerance.
    REQUIRE(e.seconds() == ref);
    REQUIRE_THAT(e.fraction().val(), WithinAbs(0.0, 1.0e-15));
  }

  SECTION("sub-second input is carried in the fraction, not rounded into the epoch") {
    Epoch e = Epoch::FromGregorian(2020, 1, 1, 12, 34, Real(56.789), Time::TDB);
    REQUIRE(e.seconds() == ExactSeconds(2020, 1, 1, 12, 34, 56));
    REQUIRE_THAT(e.fraction().val(), WithinAbs(0.789, 1.0e-12));
  }

  SECTION("the string form agrees with the numeric form") {
    Epoch a = Epoch::FromGregorian("2020-07-04T12:34:56", Time::TDB);
    Epoch b = Epoch::FromGregorian(2020, 7, 4, 12, 34, Real(56.0), Time::TDB);
    REQUIRE(a.seconds() == b.seconds());
    REQUIRE_THAT((a - b).val(), WithinAbs(0.0, 1.0e-15));
  }

  SECTION("differences between calendar dates are exact") {
    // One year apart to the second; the difference must be exactly the day count.
    Epoch a = Epoch::FromGregorian(2024, 3, 1, 7, 8, Real(9.0), Time::TAI);
    Epoch b = Epoch::FromGregorian(2025, 3, 1, 7, 8, Real(9.0), Time::TAI);
    REQUIRE_THAT((b - a).val(), WithinAbs(365.0 * 86400.0, 1.0e-9));
  }
}
