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

namespace lupnt {
namespace jpl_de {
struct EphemerisData {
  int KSIZE;
  int NCOEFF;
  std::string info;
  std::string start_info;
  std::string final_info;
  double jd_start_tbd;
  double jd_end_tbd;
  double step;
  std::vector<std::string> constant_names;
  std::vector<double> constant_values;
  std::vector<int> coeff_start_locations;
  std::vector<int> num_coefficients;
  std::vector<int> num_properties;
};

void ParseGroup1010(std::ifstream& infile, EphemerisData& data) {
  std::string line, empty_line;
  std::getline(infile, empty_line);
  std::getline(infile, data.info);
  std::getline(infile, data.start_info);
  std::getline(infile, data.final_info);
  std::getline(infile, empty_line);
}

void ParseGroup1030(std::ifstream& infile, EphemerisData& data) {
  std::string line, empty_line;
  std::getline(infile, empty_line);
  std::getline(infile, line);
  std::istringstream iss(line);
  iss >> data.jd_start_tbd >> data.jd_end_tbd >> data.step;
  std::getline(infile, empty_line);
}

void ParseGroup1040(std::ifstream& infile, EphemerisData& data) {
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

void ParseGroup1041(std::ifstream& infile, EphemerisData& data) {
  std::string line, empty_line;
  std::getline(infile, empty_line);

  int num_values;
  infile >> num_values;

  std::string value_str;
  for (int i = 0; i < num_values; ++i) {
    infile >> value_str;
    for (char& c : value_str)
      if (c == 'D') c = 'e';
    data.constant_values.push_back(std::stod(value_str));
  }
  std::getline(infile, empty_line);
  std::getline(infile, empty_line);
}

void ParseGroup1050(std::ifstream& infile, EphemerisData& data) {
  std::string line, empty_line;
  std::getline(infile, empty_line);
  std::vector<int> start_locations, num_coefficients, num_properties;
  // Read the next three lines and store the values in respective vectors
  for (int i = 0; i < 3; ++i) {
    std::getline(infile, line);
    std::istringstream iss(line);
    int value;
    while (iss >> value) {
      if (i == 0) {
        start_locations.push_back(value);
      } else if (i == 1) {
        num_coefficients.push_back(value);
      } else if (i == 2) {
        num_properties.push_back(value);
      }
    }
  }
  std::getline(infile, empty_line);

  // Assign the parsed values to the struct
  data.coeff_start_locations = start_locations;
  data.num_coefficients = num_coefficients;
  data.num_properties = num_properties;
}

void ReadEphemerisFile(const std::string& filepath, EphemerisData& data) {
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
}  // namespace jpl_de
}  // namespace lupnt