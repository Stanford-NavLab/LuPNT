#include <Eigen/Dense>
#include <chrono>
#include <iostream>

using real = double;
using Vec3 = Eigen::Vector3d;

// Function to convert geocentric coordinates to cartesian coordinates
// (pass-by-value)
Vec3 GeocentricToCartesianValue(Vec3 r_geo, real radius) {
  real r = radius;
  real theta = r_geo(1);  // Latitude
  real phi = r_geo(2);    // Longitude

  real x = r * std::cos(theta) * std::cos(phi);
  real y = r * std::cos(theta) * std::sin(phi);
  real z = r * std::sin(theta);

  return Vec3(x, y, z);
}

// Function to convert geocentric coordinates to cartesian coordinates
// (pass-by-const-reference)
Vec3 GeocentricToCartesianRef(const Vec3 &r_geo, real radius) {
  real r = radius;
  real theta = r_geo(1);  // Latitude
  real phi = r_geo(2);    // Longitude

  real x = r * std::cos(theta) * std::cos(phi);
  real y = r * std::cos(theta) * std::sin(phi);
  real z = r * std::sin(theta);

  return Vec3(x, y, z);
}

int main() {
  Vec3 r_geo(6371.0, M_PI / 4,
             M_PI / 4);  // Example geocentric coordinates (radius, theta, phi)
  real radius = 6371.0;  // Earth's radius in kilometers
  const int iterations = 100000000;  // Number of iterations for timing

  // Timing the function with pass-by-value
  auto start_value = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i) {
    Vec3 r_cartesian = GeocentricToCartesianValue(r_geo, radius);
  }
  auto end_value = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration_value = end_value - start_value;
  std::cout << "Time taken by GeocentricToCartesianValue: "
            << duration_value.count() << " seconds" << std::endl;

  // Timing the function with pass-by-const-reference
  auto start_ref = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i) {
    Vec3 r_cartesian = GeocentricToCartesianRef(r_geo, radius);
  }
  auto end_ref = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration_ref = end_ref - start_ref;
  std::cout << "Time taken by GeocentricToCartesianRef: "
            << duration_ref.count() << " seconds" << std::endl;

  return 0;
}