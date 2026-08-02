#pragma once

#include <filesystem>
#include <functional>

#include "lupnt/core/definitions.h"

namespace lupnt {

  /// Which IERS product the loaded EOP table came from.
  enum class EopSource {
    C04,     ///< IERS EOP 14 C04: final series only, retrospective, no predictions.
    Finals,  ///< IERS finals.all ("Bulletin A"): rapid solution plus ~1 year of predictions.
  };

  /// What the celestial-pole columns of the source file mean. The two finals variants share an
  /// identical fixed-width layout but not the meaning of those columns, so this cannot be
  /// inferred from the layout and must be stated.
  enum class EopNutation {
    Iau1980,  ///< dPsi, dEps -- nutation offsets in longitude/obliquity (C04, finals.all).
    Iau2000,  ///< dX, dY -- CIP offsets to the IAU 2006/2000A model (finals2000A.all).
  };

  // Struct to hold EOP data
  struct EopFileData {
    VecXi years;
    VecXi months;
    VecXi days;
    VecXd mjds_utc;
    VecXd x;            ///< Polar motion x ["]
    VecXd y;            ///< Polar motion y ["]
    VecXd ut1_utc;      ///< UT1 - UTC [s]
    VecXd lod;          ///< Length of day excess [s]
    VecXd dpsi;         ///< Nutation correction in longitude ["]
    VecXd deps;         ///< Nutation correction in obliquity ["]
    VecXd xErr;         ///< Formal error of x ["]
    VecXd yErr;         ///< Formal error of y ["]
    VecXd ut1_utc_err;  ///< Formal error of ut1_utc [s]
    VecXd lod_err;      ///< Formal error of lod [s]
    VecXd dpsi_err;     ///< Formal error of dpsi ["]
    VecXd deps_err;     ///< Formal error of deps ["]
    VecXd dX;           ///< CIP offset dX ["], IAU2000 sources only; all-zero otherwise
    VecXd dY;           ///< CIP offset dY ["], IAU2000 sources only; all-zero otherwise
    VecXd dX_err;       ///< Formal error of dX ["]
    VecXd dY_err;       ///< Formal error of dY ["]
    /// 1 where the row is not a final/measured value: a `P`-flagged Bulletin A prediction for
    /// EopSource::Finals, or a row carrying the 0.999 formal-error sentinel (IERS's marker for
    /// values that are not final) for EopSource::C04. 0 elsewhere.
    VecXi is_prediction;
    EopSource source{EopSource::C04};
    /// Which pair of columns the file carried. dpsi/deps are populated for Iau1980 and dX/dY for
    /// Iau2000; the unused pair is zero-filled rather than left empty, so indexing is always safe.
    EopNutation nutation{EopNutation::Iau1980};
  };

  /// Earth orientation parameters interpolated to a single epoch, together with the IERS formal
  /// (1-sigma) errors of the measured series.
  ///
  /// Note the sigma_* fields are the *formal errors* of the C04 solution, not celestial-pole
  /// offsets (dX/dY) and not TAI-UTC -- they were previously misnamed dx_pole/dy_pole/tai_utc.
  /// If CIP offsets are ever needed they get their own fields; TAI-UTC comes from
  /// interfaces/tai_utc.h.
  struct EopData {
    Real x_pole;         ///< Polar motion x [rad]
    Real y_pole;         ///< Polar motion y [rad]
    Real ut1_utc;        ///< UT1 - UTC [s]
    Real lod;            ///< Length of day excess [s]
    Real dpsi;           ///< Nutation correction in longitude (IAU1980) [rad]
    Real deps;           ///< Nutation correction in obliquity (IAU1980) [rad]
    Real sigma_x_pole;   ///< Formal error of x_pole [rad]
    Real sigma_y_pole;   ///< Formal error of y_pole [rad]
    Real sigma_ut1_utc;  ///< Formal error of ut1_utc [s]
    /// IERS celestial-pole offsets: observed corrections to the IAU 2006/2000A CIP, applied by
    /// RotPrecessionNutation as X += dX, Y += dY (IERS Conventions 2010 eq. 5.14). Zero unless an
    /// EopNutation::Iau2000 source is loaded. They are ~0.3 mas rms (~1 cm at the Earth's
    /// surface), so they matter for cm-level work and not for anything coarser.
    Real dX;  ///< [rad]
    Real dY;  ///< [rad]
  };

  /// Additive perturbation of the Earth-orientation parameters, in the sense
  /// `perturbed = nominal + delta`. Injected inside GetEopData, so it reaches *everything* that
  /// consumes EOP -- polar motion, the UT1 time conversion and hence sidereal rotation, LOD, and
  /// the celestial-pole offsets -- without a parallel frame stack.
  ///
  /// Members are `Real` rather than `double` on purpose: an autodiff-seeded perturbation keeps
  /// its derivative components when copied into the global, so d(observable)/d(EOP) falls out of
  /// an ordinary forward evaluation. Zero-initialised, so a default-constructed perturbation is
  /// exactly "no perturbation".
  struct EopPerturbation {
    Real dx_pole{0};  ///< Polar motion x [rad]
    Real dy_pole{0};  ///< Polar motion y [rad]
    Real dut1{0};     ///< UT1-UTC [s]
    Real dlod{0};     ///< Length of day [s]
    Real ddX{0};      ///< CIP offset dX [rad]
    Real ddY{0};      ///< CIP offset dY [rad]
  };

  /// @brief Apply a constant EOP perturbation to every subsequent GetEopData call.
  ///
  /// The mechanism behind the Step-2 sensitivity work and behind truth-vs-model Monte Carlo:
  /// generate truth with a perturbation applied and run the filter without it, or difference two
  /// evaluations to get a numerical partial. Replaces any previously-set perturbation or
  /// perturbation function.
  void SetEopPerturbation(const EopPerturbation& perturbation);

  /// @brief Apply a time-varying EOP perturbation, evaluated at each requested epoch.
  ///
  /// Used to inject a realistic prediction-error time series (autocorrelated, seasonal) as truth,
  /// rather than a constant bias. `fn` receives the epoch [MJD, UTC] and must be thread-safe: it
  /// is called under the EOP lock from whatever threads run frame conversions.
  void SetEopPerturbationFunction(std::function<EopPerturbation(Real mjd_utc)> fn);

  /// @brief Remove any EOP perturbation, restoring the unperturbed table.
  void ClearEopPerturbation();

  /// @brief The perturbation in force at `mjd_utc` (zero-valued when none is set).
  EopPerturbation GetEopPerturbation(Real mjd_utc);

  /// @brief Whether anything can currently make the celestial-pole offsets non-zero -- either an
  /// IAU2000 source is loaded or a perturbation is active.
  ///
  /// RotPrecessionNutation consults this to skip an EOP table lookup it would otherwise have to
  /// do on every call; the offsets are zero for the default C04 source, so the common path stays
  /// exactly as cheap as it was.
  bool EopHasCelestialPoleOffsets();

  /// Epoch range spanned by the loaded EOP table [MJD, UTC].
  struct EopCoverage {
    double mjd_first;
    double mjd_last;
    /// Last epoch backed by a final/measured value. Beyond this and up to mjd_last the table is
    /// predicted (Bulletin A `P` rows, or C04's 0.999-sentinel rows); beyond mjd_last GetEopData
    /// clamps. Equals mjd_last when the table carries no predicted rows.
    double mjd_last_measured;
    EopSource source;
  };

  /// @brief Declare which EOP product LuPNT should use, resolved when the table is first needed.
  ///
  /// The EOP table is global, lazily-loaded state: without this, the source is whatever the first
  /// caller happened to load, and a `Load*(path, force=false)` issued after anything has already
  /// touched Earth orientation (any frame conversion) is silently ignored -- so a run can believe
  /// it opted into one product while carrying another. This declares the choice instead of
  /// racing for it: call it at startup, or let it default.
  ///
  /// Ordering does not matter. If nothing is loaded yet the choice is recorded and applied on
  /// first use; if a table is already loaded and differs, it is reloaded immediately.
  ///
  /// The default is EopSource::C04, which is a deliberate reproducibility choice rather than an
  /// accuracy one: C04 is settled, so a result is repeatable, whereas the finals product is
  /// re-issued daily as predictions become rapid and rapid becomes final. C04 is *retrospective*
  /// though, lagging by months, and GetEopData holds its last row constant beyond that -- an
  /// arbitrary error (~20 m of Earth-surface displacement at the time of writing, dominated by
  /// UT1) that does not shrink just because the file is recent. Any run at a present-day or
  /// future epoch should select EopSource::Finals and check GetEopCoverage().
  ///
  /// @param source   Which product to use.
  /// @param filepath Optional explicit file. When empty, the default for the source is used:
  ///                  the bundled GetFilePath(EOP_FILENAME) for C04, or
  ///                  GetDataPath()/"planetary_coeff"/EOP_FINALS_FILENAME for Finals (written by
  ///                  LoadLatestEopFinalsFromIers -- nothing is bundled, so download or supply
  ///                  one first). Throws if the file does not exist, at the point of the call
  ///                  rather than later inside a frame conversion.
  void SetEopSource(EopSource source, const std::filesystem::path& filepath = {});

  /// @brief The currently declared EOP source (see SetEopSource), which may not be loaded yet.
  ///
  /// For what is actually loaded, use GetEopCoverage().source -- the two agree once the table has
  /// been read, and any direct Load* call updates this selection to match.
  EopSource GetEopSource();

  /// @brief Parse an IERS EOP 14 C04 (IAU1980) time-series file into the global EOP table.
  ///
  /// Called lazily (via GetEopFileData/GetEopData) the first time Earth orientation data is
  /// needed, using the bundled file at GetFilePath(EOP_FILENAME); also called by
  /// LoadLatestEopFromIers after downloading a fresh file. Subsequent calls are no-ops unless
  /// `force` is set, since all later frame conversions (GcrfToItrf/ItrfToGcrf via
  /// RotPolarMotion/RotSideralMotion) and time conversions (UT1-UTC) read from this table.
  ///
  /// @param filepath Path to an IERS EOP 14 C04 (IAU1980) text file (14 header lines followed by
  ///                  whitespace-separated columns: year, month, day, MJD (UTC), x_pole, y_pole,
  ///                  UT1-UTC, LOD, dpsi, deps, and their formal errors).
  /// @param force    If true, reload and replace any previously-loaded EOP table even if data is
  ///                  already present (default: false).
  void LoadEopFileData(const std::filesystem::path& filepath, bool force = false);

  /// @brief Download the latest IERS EOP 14 C04 (IAU1980) time series and load it, replacing any
  /// previously-loaded EOP data.
  ///
  /// Used at simulation startup (or on demand) to refresh the polar-motion/UT1-UTC table used by
  /// GetEopData/GetUt1UtcDifference with up-to-date IERS measurements instead of the bundled
  /// snapshot. Caches the downloaded file under
  /// GetDataPath()/"planetary_coeff"/"eopc04_iers_latest.txt". On download failure, ensures the
  /// bundled EOP file is loaded (if not already) and returns false.
  ///
  /// @param force If true, reload even if EOP data has already been loaded (default: true).
  /// @return      True if the latest IERS data was downloaded and loaded successfully, false if
  ///              the download failed (in which case the bundled data is used as a fallback).
  bool LoadLatestEopFromIers(bool force = true);

  /// @brief Parse an IERS finals ("Bulletin A") EOP file into the global EOP table, replacing any
  /// previously-loaded EOP data.
  ///
  /// This is the only loader that reaches the present day and beyond: C04 is a retrospective
  /// final series that lags by months, whereas finals files carry the rapid solution up to within
  /// a few days of issue plus roughly a year of `P`-flagged predictions. Use it wherever a
  /// simulation epoch is near or after the C04 cutoff -- otherwise GetEopData holds the last C04
  /// row constant, which is an arbitrary (not merely stale) error dominated by UT1.
  ///
  /// Expects the fixed-width IAU1980 variant (`finals.all` / `finals.all.iau1980.txt`), whose
  /// nutation columns are dPsi/dEpsilon, matching EopFileData's dpsi/deps. The IAU2000 variant
  /// (`finals2000A.all`) shares the layout but carries celestial-pole offsets dX/dY in those
  /// columns, so loading it fills dpsi/deps with dX/dY; polar motion, UT1-UTC and LOD are
  /// identical between the variants and nothing downstream currently consumes dpsi/deps.
  ///
  /// Finals rows drop trailing columns as data runs out -- LOD typically stops a day or so before
  /// the last polar-motion row and the nutation columns stop partway into the prediction span, and
  /// the file ends with date-only rows. Rows without polar motion and UT1-UTC are skipped; for
  /// rows missing only the optional columns, the last valid LOD/nutation value is held constant
  /// (a load-time Logger::Info reports where each ran out).
  ///
  /// @param filepath  Path to an IERS finals file.
  /// @param force     If true, reload and replace any previously-loaded EOP table even if data is
  ///                   already present (default: false).
  /// @param nutation  What the celestial-pole columns mean. The layout is identical between the
  ///                   variants, so this cannot be detected and getting it wrong silently files
  ///                   dX/dY under dpsi/deps or vice versa. Use Iau2000 for `finals2000A.all`.
  void LoadEopFinalsFileData(const std::filesystem::path& filepath, bool force = false,
                             EopNutation nutation = EopNutation::Iau1980);

  /// @brief Download the latest IERS finals ("Bulletin A", IAU1980) file and load it, replacing
  /// any previously-loaded EOP data.
  ///
  /// The counterpart to LoadLatestEopFromIers for the finals product; prefer this one for any run
  /// at a present-day or future epoch. Caches the downloaded file under
  /// GetDataPath()/"planetary_coeff"/"finals.all.iau1980.txt". On download failure, leaves any
  /// already-loaded table alone, ensures the bundled C04 file is loaded if nothing was, and
  /// returns false.
  ///
  /// @param force If true, reload even if EOP data has already been loaded (default: true).
  /// @return      True if the latest finals file was downloaded and loaded successfully.
  bool LoadLatestEopFinalsFromIers(bool force = true);

  /// @brief Interpolate Earth orientation parameters (polar motion, UT1-UTC, LOD, nutation
  /// corrections) to a requested epoch from the loaded IERS EOP table.
  ///
  /// Called by frame-conversion routines such as RotPolarMotion, RotSideralMotion, and
  /// RotSideralMotionDot (see frame_conversions.cc) to obtain the EOP values needed for the
  /// GCRF<->ITRF rotation at a given epoch. Lazily loads the bundled EOP file via
  /// LoadEopFileData/GetFilePath(EOP_FILENAME) on first use. For epochs outside the table's
  /// range, the nearest endpoint's values are held constant (no extrapolation); otherwise a
  /// 3rd-order Lagrange interpolation is used.
  ///
  /// Clamping is silent in the sense that the returned values carry no flag, but the first clamp
  /// after each table load emits a Logger::Warn naming the requested epoch and the table
  /// coverage: holding a years-stale UT1-UTC constant is a large, easily-overlooked Earth
  /// orientation error (a ~1 ms UT1 error is ~46 cm at the Earth's surface). The first epoch that
  /// falls in a predicted span (past EopCoverage::mjd_last_measured) likewise warns once. Use
  /// GetEopCoverage() to check an epoch up front.
  ///
  /// @param mjd_utc Epoch at which to evaluate the EOP [MJD, UTC]
  /// @return        Interpolated EOP values (x_pole/y_pole [rad], ut1_utc [s], lod [s],
  ///                dpsi/deps [rad], sigma_x_pole/sigma_y_pole [rad], sigma_ut1_utc [s])
  EopData GetEopData(Real mjd_utc);

  /// @brief Return the epoch range covered by the loaded EOP table, loading the bundled file
  /// first if no EOP data has been loaded yet.
  ///
  /// Outside this range GetEopData holds the nearest endpoint constant rather than
  /// extrapolating, so callers running epochs beyond the end of the file (a common case for
  /// future-epoch simulations against a snapshot table) can detect and quantify that up front.
  ///
  /// @return {mjd_first, mjd_last} of the loaded table [MJD, UTC]
  EopCoverage GetEopCoverage();

  /// @brief Return a pointer to the raw loaded IERS EOP table, loading the bundled file first if
  /// no EOP data has been loaded yet.
  ///
  /// Used where direct access to the full raw time series (e.g. for plotting or diagnostics) is
  /// needed, rather than the interpolated single-epoch values returned by GetEopData.
  ///
  /// @return Pointer to the global EopFileData table (owned internally; non-null after the call).
  EopFileData* GetEopFileData();

  /// @brief Get the UT1-UTC time difference at a given epoch.
  ///
  /// Used by time-conversion routines (see conversions/time_conversions.cc) to convert between
  /// UTC and UT1 (e.g. for Earth rotation angle / sidereal time computations). Thin wrapper
  /// around GetEopData that extracts just the ut1_utc field.
  ///
  /// @param mjd_utc Epoch at which to evaluate UT1-UTC [MJD, UTC]
  /// @return        UT1 - UTC [s]
  Real GetUt1UtcDifference(Real mjd_utc);

}  // namespace lupnt
