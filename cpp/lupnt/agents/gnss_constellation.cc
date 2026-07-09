#include "lupnt/agents/gnss_constellation.h"

#include <algorithm>
#include <cmath>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"
#include "lupnt/interfaces/antex_loader.h"
#include "lupnt/interfaces/sp3_loader.h"
#include "lupnt/measurements/comms_utils.h"
#include "lupnt/numerics/interpolation.h"

namespace lupnt {

  namespace {
    ChebyshevFitModel FitStateHistoryChebyshev(const VecXd& t_tai, const MatXd& rv_eci) {
      double t_start = t_tai(0);
      double t_end = t_tai(t_tai.size() - 1);
      double min_step = t_end - t_start;
      for (int i = 1; i < t_tai.size(); i++) {
        min_step = std::min(min_step, t_tai(i) - t_tai(i - 1));
      }

      return FitChebyshevModel(
          [&t_tai, &rv_eci](double t) {
            VecXd state(6);
            for (int k = 0; k < 6; k++) {
              VecXd col = rv_eci.col(k);
              state(k) = LinearInterp1d(t_tai, col, t);
            }
            return state;
          },
          t_start, t_end, 6, std::max(1.0, min_step), 2);
    }
  }  // namespace

  // Constructors ***************************************************************

  GnssConstellation::GnssConstellation(GnssConst gnss_const) : gnss_const_(gnss_const) {
    name_ = std::string(enum_name(gnss_const_));
  }

  GnssConstellation::GnssConstellation(Config& config) {
    gnss_const_ = enum_cast<GnssConst>(config["gnss_const"].as<std::string>()).value();
    name_ = config["name"] ? config["name"].as<std::string>() : std::string(enum_name(gnss_const_));
  }

  // Setup **********************************************************************

  void GnssConstellation::SetSatelliteStates(const std::vector<int>& prns, const VecXd& t_tai,
                                             const std::vector<MatXd>& rv_eci) {
    LUPNT_CHECK(prns.size() == rv_eci.size(), "Number of PRNs must match number of state histories",
                "GnssConstellation");
    LUPNT_CHECK(t_tai.size() >= 2, "Ephemeris must contain at least 2 epochs", "GnssConstellation");

    prns_ = prns;
    t_tai_ = t_tai;
    rv_eci_.clear();
    rv_eci_cheby_.clear();
    for (size_t i = 0; i < prns.size(); i++) {
      LUPNT_CHECK(rv_eci[i].rows() == t_tai.size() && rv_eci[i].cols() == 6,
                  "Each satellite state history must be [N_epochs x 6]", "GnssConstellation");
      rv_eci_[prns[i]] = rv_eci[i];
      rv_eci_cheby_[prns[i]] = FitStateHistoryChebyshev(t_tai_, rv_eci[i]);
    }
  }

  void GnssConstellation::LoadEphemeris(const std::filesystem::path& filepath) {
    LUPNT_CHECK(std::filesystem::exists(filepath),
                "GNSS ephemeris file not found: " + filepath.string(), "GnssConstellation");

    H5Easy::File file(filepath.string(), H5Easy::File::ReadOnly);
    LUPNT_CHECK(
        file.exist("/prns") && file.exist("/t_tai"),
        "Ephemeris file is missing required datasets '/prns' and/or '/t_tai': " + filepath.string(),
        "GnssConstellation");

    std::vector<int> prns = H5Easy::load<std::vector<int>>(file, "/prns");
    VecXd t_tai = H5Easy::load<VecXd>(file, "/t_tai");

    std::vector<MatXd> rv_eci;
    rv_eci.reserve(prns.size());
    for (int prn : prns) {
      std::string dset = "/rv_eci/" + std::to_string(prn);
      LUPNT_CHECK(file.exist(dset),
                  "Ephemeris file is missing dataset '" + dset + "': " + filepath.string(),
                  "GnssConstellation");
      rv_eci.push_back(H5Easy::load<MatXd>(file, dset));
    }

    SetSatelliteStates(prns, t_tai, rv_eci);
  }

  void GnssConstellation::SaveEphemeris(const std::filesystem::path& filepath) const {
    LUPNT_CHECK(!prns_.empty(), "No satellite states to save (call SetSatelliteStates first)",
                "GnssConstellation");

    H5Easy::File file(filepath.string(), H5Easy::File::Overwrite);
    H5Easy::dump(file, "/prns", prns_);
    H5Easy::dump(file, "/t_tai", t_tai_);
    for (int prn : prns_) {
      H5Easy::dump(file, "/rv_eci/" + std::to_string(prn), rv_eci_.at(prn));
    }
    file.flush();
  }

  void GnssConstellation::SetupSatelliteStatesFromFiles(
      const std::vector<std::filesystem::path>& sp3_filepaths,
      const std::filesystem::path& antex_filepath, const VecXd& t_tai, GnssFreq freq,
      const std::vector<int>& prns_in) {
    LUPNT_CHECK(t_tai.size() >= 2, "Ephemeris must contain at least 2 epochs", "GnssConstellation");
    LUPNT_CHECK(!sp3_filepaths.empty(), "At least one SP3 file path must be provided",
                "GnssConstellation");

    Sp3Loader sp3(sp3_filepaths);
    AntexLoader antex(antex_filepath);

    std::string letter = AntexLoader::GnssLetter(gnss_const_);

    // Determine the PRN list: explicit, or every PRN of `gnss_const_` present
    // in the loaded SP3 ephemeris (and -- to avoid runtime errors -- also
    // present in the ANTEX file).
    std::vector<int> prns = prns_in;
    if (prns.empty()) {
      for (const std::string& sat_id : sp3.GetSatellites()) {
        if (sat_id.size() == 3 && sat_id[0] == letter[0]) {
          int prn = std::stoi(sat_id.substr(1));
          if (antex.HasSatellite(gnss_const_, prn)) prns.push_back(prn);
        }
      }
      std::sort(prns.begin(), prns.end());
    }
    LUPNT_CHECK(!prns.empty(),
                "No PRNs of the requested constellation found in the SP3/ANTEX files",
                "GnssConstellation");

    int n_epochs = static_cast<int>(t_tai.size());
    std::vector<MatXd> rv_eci_list;
    rv_eci_list.reserve(prns.size());

    for (int prn : prns) {
      std::string sat_id = AntexLoader::SatId(gnss_const_, prn);
      MatXd rv_eci(n_epochs, 6);

      for (int k = 0; k < n_epochs; k++) {
        Real t = t_tai(k);

        // SP3 (center-of-mass) ECEF position/velocity
        Vec6 rv_sp3_ecef = sp3.GetPosVel(sat_id, t);
        Vec3d pos_sp3_ecef(rv_sp3_ecef(0).val(), rv_sp3_ecef(1).val(), rv_sp3_ecef(2).val());

        // Antenna phase-center offset (NEU, meters) -> apply correction in ECEF. ANTEX omits
        // frequencies a satellite does not transmit (e.g. L5 on older GPS blocks); fall back to
        // a zero offset there -- such links are dropped later by the C/N0 threshold anyway.
        Vec3d pco_neu_m = antex.HasPco(gnss_const_, prn, freq, t)
                              ? antex.GetPco(gnss_const_, prn, freq, t)
                              : Vec3d::Zero();
        Vec3d pos_corrected_ecef = AntexLoader::ApplyPcoCorrectionEcef(t, pos_sp3_ecef, pco_neu_m);

        // Retain the (uncorrected) SP3 velocity -- the PCO is constant in the
        // satellite body frame, so its ECEF-velocity contribution is negligible.
        Vec6 rv_corrected_ecef;
        rv_corrected_ecef << Real(pos_corrected_ecef(0)), Real(pos_corrected_ecef(1)),
            Real(pos_corrected_ecef(2)), rv_sp3_ecef(3), rv_sp3_ecef(4), rv_sp3_ecef(5);

        Real t_tdb = ConvertTime(t, Time::TAI, Time::TDB);
        Vec6 rv_eci_k = ConvertFrame(t_tdb, rv_corrected_ecef, Frame::ECEF, Frame::ECI, false);
        for (int c = 0; c < 6; c++) rv_eci(k, c) = rv_eci_k(c).val();
      }

      rv_eci_list.push_back(rv_eci);
    }

    SetSatelliteStates(prns, t_tai, rv_eci_list);
  }

  void GnssConstellation::SetupTransmitters() {
    LUPNT_CHECK(!prns_.empty(),
                "Satellite states must be set (SetSatelliteStates / LoadEphemeris) "
                "before SetupTransmitters",
                "GnssConstellation");

    antennas_.clear();
    P_tx_.clear();
    freq_list_.clear();

    for (int prn : prns_) {
      // Reuse the existing GnssTransmitter device, which already loads the
      // antenna gain pattern(s) and transmit power(s) for the given GNSS
      // constellation and PRN (see `GnssTransmitter::InitGps` / `InitGalileo`
      // / `InitQzss` in `cpp/lupnt/devices/space_comms.cc`).
      GnssTransmitter tx(gnss_const_, prn);
      antennas_[prn] = tx.GetAntennas();
      P_tx_[prn] = tx.GetTransmitPower();
      if (freq_list_.empty()) freq_list_ = tx.GetFreqList();
    }
  }

  // Queries ********************************************************************

  bool GnssConstellation::IsFaultPrn(int prn) const {
    return std::find(fault_prns_.begin(), fault_prns_.end(), prn) != fault_prns_.end();
  }

  Vec6 GnssConstellation::GetSatelliteStateEci(int prn, Real t_tai) const {
    auto it = rv_eci_.find(prn);
    LUPNT_CHECK(it != rv_eci_.end(), "PRN " + std::to_string(prn) + " not found in constellation",
                "GnssConstellation");

    auto it_cheby = rv_eci_cheby_.find(prn);
    LUPNT_CHECK(it_cheby != rv_eci_cheby_.end(),
                "Chebyshev ephemeris missing for PRN " + std::to_string(prn), "GnssConstellation");

    VecX state_x;
    LUPNT_CHECK(it_cheby->second.Eval(t_tai, &state_x, nullptr),
                "Ephemeris interpolation time is out of range", "GnssConstellation");
    Vec6 state = state_x;
    return state;
  }

  Real GnssConstellation::ComputeRange(int prn, const Vec3& r_rx_eci, Real t_tai) const {
    Vec6 rv_tx = GetSatelliteStateEci(prn, t_tai);
    return (r_rx_eci - rv_tx.head(3)).norm();
  }

  Real GnssConstellation::ComputeRangeRate(int prn, const Vec6& rv_rx_eci, Real t_tai) const {
    Vec6 rv_tx = GetSatelliteStateEci(prn, t_tai);
    Vec3 dr = rv_rx_eci.head(3) - rv_tx.head(3);
    Vec3 dv = rv_rx_eci.tail(3) - rv_tx.tail(3);
    Real range = dr.norm();
    return dr.dot(dv) / range;
  }

  const MatXd& GnssConstellation::GetSatelliteStateHistoryEci(int prn) const {
    auto it = rv_eci_.find(prn);
    LUPNT_CHECK(it != rv_eci_.end(), "PRN " + std::to_string(prn) + " not found in constellation",
                "GnssConstellation");
    return it->second;
  }

  const ChebyshevFitModel& GnssConstellation::GetSatelliteStateChebyshevEci(int prn) const {
    auto it = rv_eci_cheby_.find(prn);
    LUPNT_CHECK(it != rv_eci_cheby_.end(),
                "Chebyshev ephemeris missing for PRN " + std::to_string(prn), "GnssConstellation");
    return it->second;
  }

  bool GnssConstellation::HasTransmitterInfo(int prn, GnssFreq freq) const {
    auto it_ant = antennas_.find(prn);
    auto it_pow = P_tx_.find(prn);
    if (it_ant == antennas_.end() || it_pow == P_tx_.end()) return false;
    return it_ant->second.find(freq) != it_ant->second.end()
           && it_pow->second.find(freq) != it_pow->second.end();
  }

  const Antenna& GnssConstellation::GetTransmitterAntenna(int prn, GnssFreq freq) const {
    LUPNT_CHECK(HasTransmitterInfo(prn, freq),
                "Antenna metadata not available for PRN " + std::to_string(prn),
                "GnssConstellation");
    return antennas_.at(prn).at(freq);
  }

  Real GnssConstellation::GetTransmitPowerDbw(int prn, GnssFreq freq) const {
    LUPNT_CHECK(HasTransmitterInfo(prn, freq),
                "Transmit-power metadata not available for PRN " + std::to_string(prn),
                "GnssConstellation");
    return P_tx_.at(prn).at(freq);
  }

  // Measurement noise models ***************************************************

  // Pseudorange (DLL) noise: direct port of `compute_gnss_pseudorange_noise`,
  // reusing the case-based formula already implemented in `SigmaDll`
  // (comms_utils.cc; the case boundaries and per-case expressions are
  // identical -- `Tc = 1/R_c` makes `1/(Tc*B_fe) = R_c/B_fe`), then applying
  // the `* (C * Tc)` chip-to-meters scaling that the Python version applies
  // on top of the dimensionless code-tracking jitter.
  Real GnssConstellation::ComputeSigmaRange(Real cn0_dbhz, GnssFreq freq) const {
    Real CN0_w = DecibelToDecimal(cn0_dbhz);
    auto it = GNSS_RC_MAP.find(freq);
    LUPNT_CHECK(it != GNSS_RC_MAP.end(), "Unknown GNSS frequency", "GnssConstellation");
    Real Rc = it->second;
    Real Tc = 1.0 / Rc;

    DllParams params;
    params.B_dll = rx_params_.Bn;
    params.T_i = rx_params_.T;
    params.d = rx_params_.D;
    params.Tc = Tc;
    params.B_fe = rx_params_.b * Rc;

    return SigmaDll(params, CN0_w) * (Real(C) * Tc);
  }

  // Pseudorange-rate (FLL) noise: direct port of
  // `compute_gnss_pseudorangerate_noise`.
  Real GnssConstellation::ComputeSigmaRangeRate(Real cn0_dbhz, GnssFreq freq) const {
    Real CN0_w = DecibelToDecimal(cn0_dbhz);
    auto it = GNSS_FREQ_MAP.find(freq);
    LUPNT_CHECK(it != GNSS_FREQ_MAP.end(), "Unknown GNSS frequency", "GnssConstellation");
    Real lambda = Real(C) / it->second;

    Real Bf = rx_params_.Bf;
    Real T = rx_params_.T;

    return lambda / (2.0 * PI * T) * sqrt(4.0 * Bf / CN0_w * (1.0 + 1.0 / (T * CN0_w)));
  }

  // Carrier-phase (PLL) noise: direct port of
  // `compute_gnss_carrier_phase_noise`. The dimensionless phase-jitter term
  // `sqrt(Bp/CN0 * (1 + 1/(2*T*CN0)))` is identical to `SigmaPll`.
  Real GnssConstellation::ComputeSigmaCarrierPhase(Real cn0_dbhz, GnssFreq freq) const {
    Real CN0_w = DecibelToDecimal(cn0_dbhz);
    auto it = GNSS_FREQ_MAP.find(freq);
    LUPNT_CHECK(it != GNSS_FREQ_MAP.end(), "Unknown GNSS frequency", "GnssConstellation");
    Real lambda = Real(C) / it->second;

    PllParams params;
    params.B_pll = rx_params_.Bp;
    params.T_i = rx_params_.T;

    return lambda / (2.0 * PI) * SigmaPll(params, CN0_w);
  }

  // Constellation interface ****************************************************

  void GnssConstellation::Setup() {
    if (antennas_.empty() && !prns_.empty()) SetupTransmitters();
  }

  void GnssConstellation::Step(Real t) {
    (void)t;  // GNSS ephemerides are precomputed; nothing to propagate here.
  }

}  // namespace lupnt
