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
#pragma once
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

  double GetTtTdbDifference(double t_tai);
  Vec6 GetLunarMantleData(Real t_tai);
  MatX6 GetLunarMantleData(VecX t_tai);

  Vec6 GetBodyPosVel(Real t_tai, NaifId target, Frame frame);
  Vec6 GetBodyPosVel(Real t_tai, NaifId center, NaifId target, Frame frame);
  MatX6 GetBodyPosVel(const VecX& t_tai, NaifId target, Frame frame);
  MatX6 GetBodyPosVel(const VecX& t_tai, NaifId center, NaifId target, Frame frame);
}  // namespace lupnt
