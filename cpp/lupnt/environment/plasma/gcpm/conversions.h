/**
 * @file gcpm/conversions.h
 * @author Keidai Iiyama
 * @brief  This file contains the definition of the coordinate conversion
 * functions
 * @version 0.1
 * @date 2025-02-07
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <array>

namespace pecsim {

  enum Frame { GEI, GEO, GSE, GSM, SM };

  void pol_to_cart(double alatr, double along, double r, std::array<double, 3>& pos_sm);

  void cart_to_pol(const std::array<double, 3>& pos_geo, double& blatr, double& blong,
                   double& rtemp);

  void sm_to_geo(const std::array<int, 2>& itime, const std::array<double, 3>& pos_sm,
                 std::array<double, 3>& pos_geo);

  void geo_to_sm(const std::array<int, 2>& itime, const std::array<double, 3>& pos_geo,
                 std::array<double, 3>& pos_sm);

  void gei_to_sm(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
                 std::array<double, 3>& x_out);

  double to_mjd(const std::array<int, 2>& itime, double& iyr, double& iday, double& ut);

  void t1(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse);

  void t2(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse);

  void t3(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse);

  void t4(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse);

  void t5(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse);

  void rotate_x(double angle, const std::array<double, 3>& x_in, std::array<double, 3>& x_out);

  void rotate_y(double angle, const std::array<double, 3>& x_in, std::array<double, 3>& x_out);

  void rotate_z(double angle, const std::array<double, 3>& x_in, std::array<double, 3>& x_out);

  void get_q_c(const std::array<int, 2>& itime, std::array<double, 3>& q_c);
}  // namespace pecsim
