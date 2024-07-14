/*
Greg Miller (gmiller@gregmiller.net) 2022
Released as public domain
http://www.celestialprogramming.com/

Class to read binary versions of JPL's Development Ephemeris.  Files in
the propper format can be obtained from:
ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux

#    Properties       Units          Center Description
0    x,y,z            km             SSB    Mercury
1    x,y,z            km             SSB    Venus
2    x,y,z            km             SSB    Earth-Moon barycenter
3    x,y,z            km             SSB    Mars
4    x,y,z            km             SSB    Jupiter
5    x,y,z            km             SSB    Saturn
6    x,y,z            km             SSB    Uranus
7    x,y,z            km             SSB    Neptune
8    x,y,z            km             SSB    Pluto
9    x,y,z            km             Earth  Moon (geocentric)
10   x,y,z            km             SSB    Sun
11   dPsi,dEps        radians               Earth Nutations in lon and obliquity
12   phi,theta,psi    radians               Lunar mantle libration
13   Ox,Oy,Oz         radians/day           Lunar mantle angular velocity
14   t                seconds               TT-TDB (at geocenter)

Example: (prints x coordinate of venus using first JD available)

JPLDE.DE de = new JPLDE.DE(@"E:\Astronomy\_Ephemeris\JPLDEBinaries\jpleph.405");
Console.WriteLine(de.getPlanet(1, de.getHeader().jdStart)[0]);

24857048.3412405
*/

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/physics/spice_interface.h"
#include "lupnt/physics/time_converter.h"

#define MAXCOEFF 1020

namespace lupnt {
namespace jpl_de {

namespace EphemID {
size_t MERCURY = 0;
size_t VENUS = 1;
size_t EARTH_MOON_BARYCENTER = 2;
size_t MARS = 3;
size_t JUPITER = 4;
size_t SATURN = 5;
size_t URANUS = 6;
size_t NEPTUNE = 7;
size_t PLUTO = 8;
size_t MOON = 9;
size_t SUN = 10;
size_t EARTH_NUTATIONS = 11;
size_t MOON_MANTLE_LIBRATIONS = 12;
size_t MOON_MANTLE_ANGULAR_VELOCITY = 13;
size_t TT_TDB = 14;
};  // namespace EphemID

struct EphemerisHeaderData {
  int KSIZE;
  int NCOEFF;
  std::string info;
  std::string start_info;
  std::string final_info;
  double jd_tdb_start;
  double jd_tdb_end;
  double step;
  std::vector<std::string> constant_names;
  std::vector<double> constant_values;
  std::vector<int> coeff_offset;
  std::vector<int> n_coeffs;
  std::vector<int> n_subintervals;
  double AU;
  double EMRAT;
  std::vector<int> n_properties = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 1};
};

struct EphemerisBlock {
  double jd_tdb_start;
  double jd_tdb_end;
  double coeff[MAXCOEFF];
};

struct EphemerisData {
  EphemerisHeaderData header;
  double jd_tdb_start;
  double jd_tdb_end;
  std::vector<EphemerisBlock> blocks;
} ephemeris_data;

double ParseDouble(const std::string& str) {
  std::string s = str;
  for (char& c : s)
    if (c == 'D') c = 'e';
  return std::stod(s);
}

void ParseGroup1010(std::ifstream& infile, EphemerisHeaderData& data) {
  std::string line, empty_line;
  std::getline(infile, empty_line);
  std::getline(infile, data.info);
  std::getline(infile, data.start_info);
  std::getline(infile, data.final_info);
  std::getline(infile, empty_line);
}

void ParseGroup1030(std::ifstream& infile, EphemerisHeaderData& data) {
  std::string line, empty_line;
  std::getline(infile, empty_line);
  std::getline(infile, line);
  std::istringstream iss(line);
  iss >> data.jd_tdb_start >> data.jd_tdb_end >> data.step;
  std::getline(infile, empty_line);
}

void ParseGroup1040(std::ifstream& infile, EphemerisHeaderData& data) {
  std::string line, empty_line;
  std::getline(infile, empty_line);
  int num_constants;
  infile >> num_constants;
  std::string constant;
  for (int i = 0; i < num_constants; ++i) {
    infile >> constant;
    data.constant_names.push_back(constant);
  }
  std::getline(infile, empty_line);
  std::getline(infile, empty_line);
}

void ParseGroup1041(std::ifstream& infile, EphemerisHeaderData& data) {
  std::string line, empty_line;
  std::getline(infile, empty_line);
  int num_values;
  infile >> num_values;
  std::string value_str;
  for (int i = 0; i < num_values; ++i) {
    infile >> value_str;
    data.constant_values.push_back(ParseDouble(value_str));
  }
  std::getline(infile, empty_line);
  std::getline(infile, empty_line);
}

void ParseGroup1050(std::ifstream& infile, EphemerisHeaderData& data) {
  std::string line, empty_line;
  std::getline(infile, empty_line);
  std::vector<int> start_locations, n_coeffs, n_subintervals;
  // Read the next three lines and store the values in respective vectors
  for (int i = 0; i < 3; ++i) {
    std::getline(infile, line);
    std::istringstream iss(line);
    int value;
    while (iss >> value) {
      if (i == 0) {
        start_locations.push_back(value);
      } else if (i == 1) {
        n_coeffs.push_back(value);
      } else if (i == 2) {
        n_subintervals.push_back(value);
      }
    }
  }
  std::getline(infile, empty_line);
  data.coeff_offset = start_locations;
  data.n_coeffs = n_coeffs;
  data.n_subintervals = n_subintervals;
}

void ReadEphemerisHeaderFile(const std::string& filepath,
                             EphemerisHeaderData& data) {
  std::ifstream infile(filepath);
  assert(infile.is_open() && "Unable to open file");
  std::string line, empty_line;
  while (std::getline(infile, line)) {
    if (line.find("KSIZE") != std::string::npos) {
      std::istringstream iss(line);
      std::string temp;
      iss >> temp >> data.KSIZE >> temp >> data.NCOEFF;
    } else if (line.find("GROUP   1010") != std::string::npos) {
      ParseGroup1010(infile, data);
    } else if (line.find("GROUP   1030") != std::string::npos) {
      ParseGroup1030(infile, data);
    } else if (line.find("GROUP   1040") != std::string::npos) {
      ParseGroup1040(infile, data);
    } else if (line.find("GROUP   1041") != std::string::npos) {
      ParseGroup1041(infile, data);
    } else if (line.find("GROUP   1050") != std::string::npos) {
      ParseGroup1050(infile, data);
    }
  }
  infile.close();
}

void ReadEphemerisCoefficientsFile(const std::string& filepath,
                                   EphemerisHeaderData& data) {
  ephemeris_data.header = data;
  ephemeris_data.blocks.clear();

  std::ifstream infile(filepath);
  assert(infile.is_open() && "Unable to open file");
  std::string value_str, line;
  int block = 1;
  int block_in, tmp;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    iss >> block_in >> tmp;
    assert(block_in == block && "Block number mismatch");
    EphemerisBlock block_data;
    infile >> value_str;
    block_data.jd_tdb_start = ParseDouble(value_str);
    infile >> value_str;
    block_data.jd_tdb_end = ParseDouble(value_str);
    for (int i = 0; i < data.NCOEFF - 2; ++i) {
      infile >> value_str;
      block_data.coeff[i] = ParseDouble(value_str);
    }
    ephemeris_data.blocks.push_back(block_data);
    std::getline(infile, line);  // Finish the line
    block++;
  }
  infile.close();
  ephemeris_data.jd_tdb_start = ephemeris_data.blocks[0].jd_tdb_start;
  ephemeris_data.jd_tdb_end = ephemeris_data.blocks.back().jd_tdb_end;
}

std::pair<double, double> ComputePolinomial(double x, const double* scale,
                                            const double* coeff, int offset,
                                            int num) {
  double x2, w0 = 0., w1 = 0., dw0 = 0., dw1 = 0., tmp;

  x = (x - scale[0]) / scale[1];
  x2 = x * 2.;
  while (--num) {
    tmp = dw1;
    dw1 = dw0;
    dw0 = w0 * 2. + dw0 * x2 - tmp;
    tmp = w1;
    w1 = w0;
    w0 = coeff[offset + num] + (x2 * w0 - tmp);
  }
  double f = coeff[offset] + (x * w0 - w1);
  double df = (w0 + x * dw0 - dw1) / scale[1];
  return {f, df};
}

VecXd GetBodyPosVel(Real t_tai, int id) {
  double t_tdb_ = ConvertTime(t_tai, "TAI", "TDB").val();
  double t_tdb = ConvertT(t_tai, TimeSys::TAI, TimeSys::TDB).val();
  double jd_tdb = t_tdb / SECS_DAY + JD_J2000;
  double jd_tdb_ = ConvertTime(t_tai, "TAI", "JDTDB").val();

  // Block
  double Dt = ephemeris_data.header.step;
  int i = int((jd_tdb - ephemeris_data.jd_tdb_start) / Dt);
  assert(i >= 0 && i < ephemeris_data.blocks.size() &&
         "Block index out of range");  // TODO: Load proper file
  EphemerisBlock& block = ephemeris_data.blocks[i];
  EphemerisHeaderData& header = ephemeris_data.header;

  // Subinterval
  double Dt_subint = Dt / header.n_subintervals[id];
  int j = int((jd_tdb - block.jd_tdb_start) / Dt_subint);
  int n_coeff = header.n_coeffs[id] * header.n_properties[id];
  int offset = header.coeff_offset[id] + j * n_coeff - 3;
  double jd_tdb_subint = block.jd_tdb_start + j * Dt_subint;

  double scale[2] = {jd_tdb_subint + Dt_subint / 2., Dt_subint / 2.};

  VecXd rv(6);
  for (int i = 0; i < 3; i++) {
    auto [pos, vel] = ComputePolinomial(jd_tdb, scale, block.coeff, offset,
                                        header.n_coeffs[id]);
    rv[i] = pos;
    rv[i + 3] = vel / SECS_DAY;
    offset += header.n_coeffs[id];
  }
  return rv;
}

}  // namespace jpl_de
}  // namespace lupnt