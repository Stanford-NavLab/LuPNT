// Builds a GroundStation object from explicit SPICE-style station metadata.
// This is useful when defining stations from code instead of loading them from
// SPICE kernels.
#include <lupnt/lupnt.h>

using namespace lupnt;

int main() {
  // Goldstone DSS-14-like location, expressed as geodetic latitude/longitude
  // and altitude before conversion to Cartesian meters.
  BodyData earth = GetBodyData(BodyId::EARTH);
  LatLonAlt lla(Vec3(35.426456, 243.110461, 1001.39), earth.fixed_frame);
  Cart3 position = LatLonAltToCart(lla, earth.R, earth.flattening);

  spice::GroundStationSpiceData data;
  data.name = "DSS14";
  data.body_id = BodyId::EARTH;
  data.frame = "ITRF";
  data.position_m = position.head(3).cast<double>();
  data.latitude_deg = lla(0).val();
  data.longitude_deg = lla(1).val();
  data.altitude_m = lla(2).val();

  GroundStation gs(data);
  std::cout << "GS Name  : " << gs.GetName() << std::endl;
  std::cout << "Lat, Lon : " << gs.GetLatitudeDouble() << " , " << gs.GetLongitudeDouble()
            << std::endl;
  std::cout << "PosVel   : " << gs.GetState().transpose() << std::endl;
}
