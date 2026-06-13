#pragma once

#include "lupnt/agents/agent.h"
#include "lupnt/conversions/coordinate_conversions.h"
#include "lupnt/environment/body.h"
#include "lupnt/interfaces/spice.h"

namespace lupnt {

  class GroundStation : public AgentWithDynamics {
  protected:
    BodyId body_id_ = BodyId::EARTH;
    Real latitude_ = 0.0;   // [deg]
    Real longitude_ = 0.0;  // [deg]
    Real altitude_ = 0.0;   // [m]
    Frame frame_ = Frame::UNDEFINED;
    int naif_id_ = 0;
    std::string spice_frame_;

    void ConfigureFromData(const spice::GroundStationSpiceData& gs_data);
    void ConfigureFromCartesian(const std::string& name, const Vec3& position_m,
                                BodyId body_id = BodyId::EARTH);
    void ConfigureFromConfig(Config& config);

  public:
    GroundStation() = default;
    GroundStation(Config& config);
    explicit GroundStation(const spice::GroundStationSpiceData& gs_data);

    static GroundStation FromSpice(Real t_tai, const std::string& station_name,
                                   BodyId center = BodyId::EARTH,
                                   const std::string& ref_frame = "ITRF93",
                                   const std::string& ab_correction = "NONE");
    static GroundStation FromSpice(Real t_tai, int station_id, BodyId center = BodyId::EARTH,
                                   const std::string& ref_frame = "ITRF93",
                                   const std::string& ab_correction = "NONE");

    double GetLatitudeDouble() const { return GetLatitudeDegDouble(); }
    double GetLongitudeDouble() const { return GetLongitudeDegDouble(); }
    double GetAltitudeDouble() const { return GetAltitudeMDouble(); }
    double GetLatitudeDegDouble() const { return latitude_.val(); }
    double GetLongitudeDegDouble() const { return longitude_.val(); }
    double GetAltitudeMDouble() const { return altitude_.val(); }
    int GetNaifId() const { return naif_id_; }
    Frame GetFrame() const { return frame_; }
    std::string GetSpiceFrame() const { return spice_frame_; }
    VecXd GetLatLonAltDouble() const { return GetLatLonAltReal().cast<double>(); }
    Real GetLatitudeReal() const { return latitude_; }
    Real GetLongitudeReal() const { return longitude_; }
    Real GetAltitudeReal() const { return altitude_; }
    VecX GetLatLonAltReal() const {
      VecX lla(3);
      lla << latitude_, longitude_, altitude_;
      return lla;
    }
  };

}  // namespace lupnt
