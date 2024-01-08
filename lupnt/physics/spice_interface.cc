/**
 * @file SpiceInterface.cpp
 * @author Stanford NAV LAB
 * @brief  SPICE Interface functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "spice_interface.h"

#include <cspice/SpiceUsr.h>
#include <cspice/SpiceZfc.h>
#include <string.h>

#include <filesystem>
#include <fstream>
#include <iostream>

#include "lupnt/core/user_file_path.h"
#include "lupnt/physics/cheby.h"

namespace lupnt {
namespace SpiceInterface {

/**
 * @brief load the Spice kernels
 *
 */
void LoadSpiceKernel(void) {
  SpiceInt kcount;

  ktotal_c("ALL", &kcount);

  if (kcount > 0) { /* return if the kernel is already loaded */
    return;
  }

  // move to cspice directory
  std::string orig_dir = std::filesystem::current_path().string();
  std::filesystem::current_path(CSPICE_KER_DIR);

  /* Load SPICE Kernels */
  // system("cd ../data/spice_kernel");
  furnsh_c("naif0012.tls");          // leap seconds
  furnsh_c("de440.bsp");             // planetary ephemeris
  furnsh_c("pck00011.tpc");          // planetary constants
  furnsh_c("moon_assoc_pa.tf");      // assing pa as default moon orientation
  furnsh_c("moon_de440_220930.tf");  // add moon_pa

  // high fidelity lunar and earth orientation parameters
  // reference:
  // http://spiftp.esac.esa.int/workshops/2012_04_ESAC_WORKSHOP/Tutorials/27_lunar-earth_pck-fk.pdf
  furnsh_c("moon_pa_de440_200625.bpc");
  furnsh_c("earth_200101_990628_predict.bpc");  // low fidelity long history
                                                // earth EOP
  furnsh_c(
      "earth_000101_230805_230512.bpc");  // shorter histrory precise earth EOP

  std::cout << "Loaded All SPICE Kernels" << std::endl;

  // Load Chebyshev coefficients
  cheby_s = spk_extract("de440.bsp", &cheby_n);
  std::cout << "Loaded Chebyshev coefficients for " + std::to_string(cheby_n) +
                   " planets."
            << std::endl;

  if (cheby_s == NULL) {
    cheby_err(
        "could not load SPK file - Please Download the SPK file. See "
        "data/ephemeris/readme.md for instructions");
  }

  // move back to original directory
  std::filesystem::current_path(orig_dir);

  return;
}

/**
 * @brief Extracts the PCK coefficients from the kernel files
 *
 */
void ExtractPckCoeffs() {
  int handle;
  SpiceInt pck_handle;
  ConstSpiceChar *pck_file = "../data/ephemeris/moon_pa_de440_200625.bpc";
  double et = 8000.0;
  int body = 301;
  double descr[5];
  char ident[40];
  int found;
  SpiceDouble record[120];
  int rsize;
  int pdeg;
  double ra;
  double dec;
  double w;
  double lambda;
  static char bref[32];
  double eulang[6];

  // The fundamental quantities defined by PCK orientation models are actually
  // Euler angles, not matrices. These Euler angles, which we call ``RA, DEC,
  // and W,'' are related to the transformation operator returned from pxform_c
  // by the equation rotate = [ W ]   [ Pi/2 - DEC ]   [ Pi/2 + RA ]
  //              3                1               3
  // To directly retrieve these angles, use the call:

  // bodeul_( &body, &et, &ra, &dec, &w, &lambda );
  // std::cout << "body: " << body << std::endl;
  // std::cout << "et: " << et << std::endl;
  // std::cout << "ra: " << ra << std::endl;
  // std::cout << "dec: " << dec << std::endl;
  // std::cout << "w: " << w << std::endl;
  // std::cout << " " << std::endl;

  pckeul_(&body, &et, &found, bref, eulang, (ftnlen)32);
  std::cout << "found:" << found << std::endl;
  std::cout << "phi: " << eulang[0] << std::endl;
  std::cout << "delta: " << eulang[1] << std::endl;
  std::cout << " " << std::endl;

  pcklof_c(pck_file, &pck_handle);  // load the PCK file
  pcksfs_(&body, &et, &handle, descr, ident, &found, (ftnlen)40);

  std::cout << "pck handle: :" << pck_handle << std::endl;
  std::cout << "handle: :" << handle << std::endl;
  std::cout << "descr: " << &descr << std::endl;
  std::cout << "ident: " << &ident << std::endl;
  std::cout << "found:" << found << std::endl;

  if (found) {
    pckr02_(&handle, descr, &et, record);
    rsize = record[1];
    pdeg = (rsize - 2) / 3 - 1;
    std::cout << "Polynomial Size:" << rsize << std::endl;
    std::cout << "Polynomial Degree:" << pdeg << std::endl;
  }

  // extract coefficients from the CSPICE PCK file
  //    pcksfs_(body, et, &handle, descr, ident, found, (ftnlen)40);
  //    pckr02_c(handle, target)
}

Vector3d GetBodyPos(std::string targetName, real epoch, std::string refFrame,
                    std::string obsName, std::string abCorrection) {
  LoadSpiceKernel();

  SpiceDouble ptarg[3];
  Vector3d targetPos;

  SpiceDouble et =
      (SpiceDouble)
          epoch.val();  // ToDO: this cuts the relatonship between et and matrix
  const char *targ =
      strcpy(new char[targetName.length() + 1], targetName.c_str());
  const char *ref = strcpy(new char[refFrame.length() + 1], refFrame.c_str());
  const char *abcorr =
      strcpy(new char[abCorrection.length() + 1], abCorrection.c_str());
  const char *obs = strcpy(new char[obsName.length() + 1], obsName.c_str());
  SpiceDouble lt;

  // void spkpos_c(ConstSpiceChar * targ, SpiceDouble et, ConstSpiceChar * ref,
  // ConstSpiceChar * abcorr,
  //               ConstSpiceChar * obs, SpiceDouble ptarg[3], SpiceDouble * lt)
  spkpos_c(targ, et, ref, abcorr, obs, ptarg, &lt);

  for (int i = 0; i < 3; i++) {
    targetPos(i) = ptarg[i];
  }

  return targetPos;
}

/**
 * @brief Get the Frame Conversion Matrix object
 *
 * @param et
 * @param from_frame
 * @param to_frame
 * @return MatrixXd
 */
Matrix6d GetFrameConversionMatrix(real et, std::string from_frame,
                                  std::string to_frame) {
  LoadSpiceKernel();
  SpiceInt bodyname;
  SpiceDouble et_spice =
      (SpiceDouble)
          et.val();  // ToDO: this cuts the relatonship between et and matrix
  double xform[6][6];
  Matrix6d M_rot;

  const char *from_frame_char =
      strcpy(new char[from_frame.length() + 1], from_frame.c_str());
  const char *to_frame_char =
      strcpy(new char[to_frame.length() + 1], to_frame.c_str());

  sxform_c(from_frame_char, to_frame_char, et_spice, xform);

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      M_rot(i, j) = xform[i][j];
    }
  }

  return M_rot;
}

/**
 * @brief Convert a string to ephemeris time
 *
 * @param str           string to be converted
 * SO (T) Formats.

   String                        Year Mon  DOY DOM  HR Min Sec
   ----------------------------  ---- ---  --- ---  -- --- ------
   1996-12-18T12:28:28           1996 Dec   na  18  12  28 28
   1986-01-18T12                 1986 Jan   na  18  12  00 00
   1986-01-18T12:19              1986 Jan   na  18  12  19 00
   1986-01-18T12:19:52.18        1986 Jan   na  18  12  19 52.18
   1986-01-18T12:19:52.18Z       1986 Jan   na  18  12  19 52.18
   1995-08T18:28:12              1995  na  008  na  18  28 12
   1995-08T18:28:12Z             1995  na  008  na  18  28 12
   1995-18T                      1995  na  018  na  00  00 00
   0000-01-01T                   1 BC Jan   na  01  00  00 00
Calendar Formats.
   String                        Year   Mon DOM  HR Min  Sec
   ----------------------------  ----   --- ---  -- ---  ------
   Tue Aug  6 11:10:57  1996     1996   Aug  06  11  10  57
   1 DEC 1997 12:28:29.192       1997   Dec  01  12  28  29.192
   2/3/1996 17:18:12.002         1996   Feb  03  17  18  12.002
   Mar 2 12:18:17.287 1993       1993   Mar  02  12  18  17.287
   1992 11:18:28  3 Jul          1992   Jul  03  11  18  28
   June 12, 1989 01:21           1989   Jun  12  01  21  00
   1978/3/12 23:28:59.29         1978   Mar  12  23  28  59.29
   17JUN1982 18:28:28            1982   Jun  17  18  28  28
   13:28:28.128 1992 27 Jun      1992   Jun  27  13  28  28.128
   1972 27 jun 12:29             1972   Jun  27  12  29  00
   '93 Jan 23 12:29:47.289       1993*  Jan  23  12  29  47.289
   27 Jan 3, 19:12:28.182        2027*  Jan  03  19  12  28.182
   23 A.D. APR 4, 18:28:29.29    0023** Apr  04  18  28  29.29
   18 B.C. Jun 3, 12:29:28.291   -017** Jun  03  12  29  28.291
   29 Jun  30 12:29:29.298       2029+  Jun  30  12  29  29.298
   29 Jun '30 12:29:29.298       2030*  Jun  29  12  29  29.298
Day of Year Formats.
   String                        Year  DOY HR Min Sec
   ----------------------------  ----  --- -- --- ------
   1997-162::12:18:28.827        1997  162 12  18 28.827
   162-1996/12:28:28.287         1996  162 12  28 28.287
   1993-321/12:28:28.287         1993  231 12  28 28.287
   1992 183// 12:18:19           1992  183 12  18 19
   17:28:01.287 1992-272//       1992  272 17  28 01.287
   17:28:01.282 272-1994//       1994  272 17  28 01.282
   '92-271/ 12:28:30.291         1992* 271 12  28 30.291
   92-182/ 18:28:28.281          1992* 182 18  28 28.281
   182-92/ 12:29:29.192          0182+ 092 12  29 29.192
   182-'92/ 12:28:29.182         1992  182 12  28 29.182
Julian Date Strings.
   jd 28272.291                  Julian Date   28272.291
   2451515.2981 (JD)             Julian Date 2451515.2981
   2451515.2981 JD               Julian Date 2451515.2981

 * @return real     ephemeris time (TDB) (seconds past the J2000 epoch)
 */
real StringToTDB(std::string str) {
  LoadSpiceKernel();
  SpiceDouble et;
  str2et_c(str.c_str(), &et);
  return et;
}

/**
 * @brief Convert string to TAI
 *
 * @param str time string
 * @return real
 */
real StringToTAI(std::string str) {
  LoadSpiceKernel();
  real et = StringToTDB(str);
  real tai = ConvertTime(et, "TDB", "TAI");
  return tai;
}

/**
 * @brief Convert string to UTC
 *
 * @param tdb time in TDB
 * @param prec precision of the output string (default 3)
 * @return std::string
 */
std::string TDBtoStringUTC(real tdb, int prec = 3) {
  LoadSpiceKernel();
  SpiceDouble et = tdb.val();
  SpiceChar str[100];
  et2utc_c(et, "C", prec, 100, str);

  return std::string(str);
}

/**
 * @brief Convert TAI to string UTC
 *
 * @param tdb time in TAI
 * @param prec precision of the output string (default 3)
 * @return std::string
 */
std::string TAItoStringUTC(real tai, int prec = 3) {
  LoadSpiceKernel();
  real et_tdb = ConvertTime(tai, "TAI", "TDB");
  std::string str = TDBtoStringUTC(et_tdb, prec);
  return str;
}

/**
 * @brief Convert time from one time system to another
 *
 * @param t        in time in seconds
 * @param from_time_type  from time system
 *  String ID   Time system
 *  ---------   --------------------------
   TAI         International Atomic Time
   TDB         Barycentric Dynamical Time
   TT          Terrestrial Time
   TDT         Terrestrial Dynamical Time (TT)
   ET          Ephemeris time, alias for TDB
   JDTDB       Julian Date relative to TDB
   JDTDT       Julian Date relative to TDT (TT)
   JED         Julian Ephemeris date (synonym to JDTDB)
   GPS         Global Positioning System Time

 * @param to_time_type  to time system
 * @return real     out time in seconds
 */
real ConvertTime(real t, std::string from_time_type, std::string to_time_type) {
  LoadSpiceKernel();
  SpiceDouble t_in = t.val();

  SpiceDouble t_out_spice =
      unitim_c(t_in, from_time_type.c_str(), to_time_type.c_str());

  double offset = t_out_spice - t_in;  // offset in seconds

  real t_out = t + offset;  // this is to convert to real

  return t_out;
}

/**
 * @brief Get the Body Position and Velocity using Chebyshev polynomials
 *
 * @param tai_MJD TAI in MJD (with the origin as JD_NOV_17_1858)
 * @param center  center body id
 * @param target  target body id
 * @return VectorX  6x1 vector of position and velocity of target body
 * in center body J2000 frame
 */
Vector6 GetBodyPosVel(const real tai, int center, int target) {
  LoadSpiceKernel();

  bool load_cent = false;
  bool load_targ = false;
  Vector6 rv_center;
  Vector6 rv_target;

  // convert TAI to TDB past J2000
  real t_s = ConvertTime(tai, "TAI", "TDB");

  // find the segment for the center and target body
  if (center == 0) {
    rv_center = Vector6::Zero();
  }
  if (target == 0) {
    rv_target = Vector6::Zero();
  }

  for (int i = 0; i < cheby_n; i++) {
    if (cheby_s[i].target == target) {
      rv_target = cheby_posvel_ad(t_s, cheby_s[i].seg, cheby_s[i].len);
      load_cent = true;
    } else if (cheby_s[i].target == center) {
      rv_center = cheby_posvel_ad(t_s, cheby_s[i].seg, cheby_s[i].len);
      load_targ = true;
    }

    if (load_cent && load_targ) {
      break;
    }
  }
  assert(load_cent && load_targ && "Chebyshev coefficients not found");

  Vector6 retState = rv_target - rv_center;
  return retState;
}

}  // namespace SpiceInterface
}  // namespace lupnt