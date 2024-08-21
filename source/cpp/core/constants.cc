#include "lupnt/core/constants.h"

namespace lupnt {
  double GetBodyRadius(NaifId body) {
    switch (body) {
      case NaifId::MOON:
        return R_MOON;
      case NaifId::EARTH:
        return R_EARTH;
      default:
        break;
    }
    throw std::runtime_error("Body radius not defined");
    return -1;
  }

  const std::ostream& operator<<(std::ostream& os, NaifId id) {
    switch (id) {
      case NaifId::SOLAR_SYSTEM_BARYCENTER:
        os << "SOLAR_SYSTEM_BARYCENTER";
        break;
      case NaifId::MERCURY_BARYCENTER:
        os << "MERCURY_BARYCENTER";
        break;
      case NaifId::VENUS_BARYCENTER:
        os << "VENUS_BARYCENTER";
        break;
      case NaifId::EARTH_MOON_BARYCENTER:
        os << "EARTH_MOON_BARYCENTER";
        break;
      case NaifId::MARS_BARYCENTER:
        os << "MARS_BARYCENTER";
        break;
      case NaifId::JUPITER_BARYCENTER:
        os << "JUPITER_BARYCENTER";
        break;
      case NaifId::SATURN_BARYCENTER:
        os << "SATURN_BARYCENTER";
        break;
      case NaifId::URANUS_BARYCENTER:
        os << "URANUS_BARYCENTER";
        break;
      case NaifId::NEPTUNE_BARYCENTER:
        os << "NEPTUNE_BARYCENTER";
        break;
      case NaifId::PLUTO_BARYCENTER:
        os << "PLUTO_BARYCENTER";
        break;
      case NaifId::SUN:
        os << "SUN";
        break;
      case NaifId::MERCURY:
        os << "MERCURY";
        break;
      case NaifId::VENUS:
        os << "VENUS";
        break;
      case NaifId::EARTH:
        os << "EARTH";
        break;
      case NaifId::MOON:
        os << "MOON";
        break;
      case NaifId::MARS:
        os << "MARS";
        break;
      case NaifId::PHOBOS:
        os << "PHOBOS";
        break;
      case NaifId::DEIMOS:
        os << "DEIMOS";
        break;
      case NaifId::JUPITER:
        os << "JUPITER";
        break;
      case NaifId::SATURN:
        os << "SATURN";
        break;
      case NaifId::URANUS:
        os << "URANUS";
        break;
      case NaifId::NEPTUNE:
        os << "NEPTUNE";
        break;
      default:
        throw std::runtime_error("Unknown NaifId");
        break;
    }

    return os;
  }

}  // namespace lupnt
