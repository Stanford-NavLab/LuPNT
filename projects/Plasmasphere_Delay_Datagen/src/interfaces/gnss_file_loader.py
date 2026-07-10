from scipy.interpolate import CubicSpline, BarycentricInterpolator
import numpy as np
from datetime import datetime, timedelta
import pylupnt as pnt
import requests
import os
import re
import gzip
from pathlib import Path
from .gnss_utils import (
    datetime_to_gpsweeks,
    datetime_to_tai,
    tai_to_datetime,
    tai_to_gps_weeks,
)


def _earthdata_session():
    session = requests.Session()
    session.trust_env = True
    username = os.environ.get("EARTHDATA_USERNAME")
    password = os.environ.get("EARTHDATA_PASSWORD")
    if username and password:
        session.auth = (username, password)
    return session


def _looks_like_html(content):
    head = content[:256].lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html")


def _download_cddis_file(fileurl, save_path):
    session = _earthdata_session()
    response = session.get(fileurl, allow_redirects=True, timeout=120)
    if response.status_code != 200:
        raise FileNotFoundError(
            f"CDDIS file request failed with HTTP {response.status_code}. URL: {fileurl}"
        )
    if _looks_like_html(response.content):
        raise PermissionError(
            "CDDIS returned an Earthdata Login HTML page instead of the requested file. "
            "Configure Earthdata credentials with ~/.netrc or EARTHDATA_USERNAME/"
            "EARTHDATA_PASSWORD, then retry. URL: {}".format(fileurl)
        )

    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    save_path.write_bytes(response.content)


class SP3Loader:
    """
    Class to load and parse precise ephemeris files in SP3 format.
    """

    def __init__(self, filenames=[], target_dt=None, sim_t=0, dt_timesys=pnt.UTC):
        """
        Initialize the SP3Loader with the filename of the SP3 file.
        """

        self.datapath = os.path.join(pnt.get_output_dir("gnss_files"), "sp3")
        if not os.path.exists(self.datapath):
            os.makedirs(self.datapath)

        is_list_filenames = isinstance(filenames, list)

        if not is_list_filenames:
            if filenames is None:
                filenames = []
            if isinstance(filenames, str):
                filenames = [filenames]

        if sim_t == 0:
            target_dts = [target_dt] if target_dt is not None else []
        else:
            target_dt_end = (
                target_dt + timedelta(seconds=sim_t) + timedelta(days=1)
            )  # margin of 1 day
            target_dt_start = target_dt - timedelta(days=1)  # margin of 1 day

            target_dt = target_dt_start
            target_dts = []
            while target_dt <= target_dt_end:
                target_dts.append(target_dt)
                target_dt += timedelta(days=1)

        if len(filenames) == 0 and len(target_dts) == 0:
            # error if no filename or date is provided
            raise ValueError("Either filename or date must be provided.")
        elif len(filenames) == 0 and len(target_dts) > 0:
            # construct filename from date
            self.filenames = [self.load_sp3(dt, dt_timesys) for dt in target_dts]
        elif len(filenames) > 0:
            self.filenames = filenames
        else:
            raise ValueError("Invalid Input")

        # store the epochs and satellites
        self.epochs = []
        self.epochs_dt = []
        self.sats = []  # use a set to avoid duplicates
        self.positions = {}  # dictionary to store positions for each satellite
        # satellites that have at least one sentinel (999999.999999 µs) clock epoch
        self.fault_clock_sats: set = set()

        # parse the SP3 file
        for filename in self.filenames:
            self.parse_sp3(filename)

        # set internal integrator
        dyn_prop = pnt.CartesianTwoBodyDynamics(pnt.GM_EARTH)
        dyn_prop.set_integrator(pnt.IntegratorType.RKF45)
        dyn_prop.set_integrator_params(
            pnt.IntegratorParams(max_iter=20, abstol=1e-10, reltol=1e-10)
        )
        self.dyn_prop = dyn_prop

    def load_sp3(self, target_dt, dt_timesys=pnt.UTC):
        # load the SP3 file from the NASA CDDIS archive
        # Construct the filename based on the date

        gps_week, day_of_week, sec_of_week, gps_dt = datetime_to_gpsweeks(target_dt, dt_timesys)

        year = gps_dt.year
        doy = gps_dt.timetuple().tm_yday

        YYYY = str(year)
        DDD = f"{doy:03d}"  # day of the year, zero-padded to 3 digits
        WWWW = f"{gps_week:04d}"  # GPS week, zero-padded to 4 digits
        DOW = f"{day_of_week:01d}"  # Day of week, zero-padded to 1 digit

        # compute gps week
        startdir = "https://cddis.nasa.gov/archive/gnss/products"

        if gps_week >= 1962:
            # COD changed the name to new format from week 1962 (2017-08-13 (day 225))
            # https://igs.org/products/#orbits_clocks
            # https://files.igs.org/pub/resource/guidelines/Guidelines_for_Long_Product_Filenames_in_the_IGS.pdf?_gl=1*1mb4zo5*_ga*ODcwMDE5MTc1LjE3NTIwMDM0MDM.*_ga_Z5RH7R682C*czE3NTIwNjg0MTAkbzgkZzEkdDE3NTIwNjg1NjYkajYwJGwwJGgw&_ga=2.140882164.791326490.1752003404-870019175.1752003403
            AAA = "COD"  # analysis center
            PPP = "MGX"  # product type
            TTT = "FIN"  # product solution type
            LEN = "01D"  # length of the product
            SMP = "05M"  # sampling resolution
            filename = f"{AAA}0{PPP}{TTT}_{YYYY}{DDD}0000_{LEN}_{SMP}_ORB.SP3"
            zfilename = f"{filename}.gz"
            fileurl = f"{startdir}/{WWWW}/{zfilename}"
            filename_unzipped = os.path.join(self.datapath, filename)
        else:
            filename = f"cod{WWWW}{DOW}.eph"
            zfilename = f"{filename}.Z"
            fileurl = f"{startdir}/{zfilename}"
            filename_unzipped = os.path.join(self.datapath, filename)

        if os.path.exists(filename_unzipped):
            print(
                f"SP3 file for GPS TIME: {gps_dt} already exists at {filename_unzipped}. Skipping download."
            )
            return filename_unzipped

        # Download the file from CDDIS
        sp3_filename = os.path.join(self.datapath, zfilename)
        _download_cddis_file(fileurl, sp3_filename)

        print(
            f"GPS Datetime: {gps_dt} | SP3 file downloaded from {fileurl} and saved as {filename_unzipped}"
        )

        # Unzip if necessary
        if sp3_filename.endswith(".gz"):
            with gzip.open(sp3_filename, "rb") as f_in:
                with open(filename_unzipped, "wb") as f_out:
                    f_out.write(f_in.read())
            os.remove(sp3_filename)
        elif sp3_filename.endswith(".Z"):
            os.system(f"uncompress {sp3_filename}")
            # Remove .Z extension if it exists
            os.remove(sp3_filename)

        print(f"SP3 file downloaded and saved as {filename_unzipped}")
        return filename_unzipped

    def parse_sp3(self, filename):
        """
        Parse SP3 file and return the satellite positions and velocities.
        """

        # positions = {}
        # epochs = []
        # epochs_dt = []
        # satellites = set()

        with open(filename, "r") as f:
            lines = f.readlines()

        current_epoch = None

        for line in lines:
            if line.startswith("*"):
                # Epoch line
                # Format: * yyyy mm dd hh mm ss.sssss
                parts = line.strip().split()
                year = int(parts[1])
                month = int(parts[2])
                day = int(parts[3])
                hour = int(parts[4])
                minute = int(parts[5])
                second = float(parts[6])

                current_epoch = pnt.convert_time(
                    pnt.gregorian_to_time(year, month, day, hour, minute, second),
                    pnt.Time.GPS,
                    pnt.TAI,
                )
                current_dt = datetime(
                    year, month, day, hour, minute, int(second), int((second - int(second)) * 1e6)
                )
                self.epochs.append(current_epoch)
                self.epochs_dt.append(current_dt)

            elif line.startswith("P"):
                # Position line
                # Format: P sv x y z clock
                sv = line[1:4].strip()
                x = float(line[4:18].strip()) * 1000  # km to meters
                y = float(line[18:32].strip()) * 1000
                z = float(line[32:46].strip()) * 1000
                clock_us = float(line[46:60].strip())  # microseconds
                # Sentinel 999999.999999 µs means clock unavailable → NaN
                clock = np.nan if clock_us > 999990.0 else clock_us * 1e-6
                if np.isnan(clock):
                    self.fault_clock_sats.add(sv)

                self.sats.append(sv)  # add satellite to the list
                self.positions.setdefault(sv, []).append((current_epoch, (x, y, z, clock)))

        self.sats = np.unique(np.array(self.sats)).tolist()  # convert to list and remove duplicates

    def get_posvelclock(self, sys, prn, epoch, out_frame=pnt.ECI, propagate=False):
        """
        Get the position of a satellite at a specific epoch.

        Parameters:
        -------------
        sat (str): Satellite identifier (e.g., 'G01' for GPS satellite 1).
        epoch (float or ndarray): Epoch in TAI seconds or a vector of epochs.
        out_frame (str): Output frame ('ECI' or 'ECEF').
        propagate (bool): Whether to propagate the orbit (set to True if the epoch is outside the range of the ephemeris data).
        """

        sat = f"{sys}{prn:02d}"  # format satellite identifier

        if isinstance(epoch, (list, np.ndarray)):
            epoch_vec = True
            if isinstance(epoch, list):
                epoch = np.array(epoch)
        else:
            epoch_vec = False

        # from scipy.interpolate import CubicSpline
        pos_sat = self.positions[sat]

        epochs_ref = np.array([e for e, _ in pos_sat])  # N x 1 array of epochs
        positions = np.array([pos for _, pos in pos_sat])  # N x 4 array (x, y, z, clock)

        # sort epochs and positions
        sorted_indices = np.argsort(epochs_ref)
        epochs_ref = epochs_ref[sorted_indices]
        positions = positions[sorted_indices]

        # Remove duplicate epochs, preferring the record with a valid (non-NaN) clock.
        # When the same epoch appears at the end of one SP3 daily file and the start of
        # the next, the first occurrence carries sentinel clock (999999.999999 µs → NaN)
        # while the second has the real value.
        keep_mask = np.ones(len(epochs_ref), dtype=bool)
        prev_idx = 0
        for i in range(1, len(epochs_ref)):
            if epochs_ref[i] == epochs_ref[prev_idx]:
                prev_nan = np.isnan(positions[prev_idx, 3])
                curr_nan = np.isnan(positions[i, 3])
                if not prev_nan and curr_nan:
                    keep_mask[i] = False  # keep previous (valid clock)
                else:
                    keep_mask[prev_idx] = False  # keep current (later record or both NaN)
                    prev_idx = i
            else:
                prev_idx = i
        epochs_ref = epochs_ref[keep_mask]
        positions = positions[keep_mask]

        # fit polynomials
        interp_type = "cubic"  # use cubic spline interpolation

        # Clock: linear interpolation on valid (non-sentinel) epochs only
        valid_clk = ~np.isnan(positions[:, 3])
        if interp_type == "cubic":
            interp_x = CubicSpline(epochs_ref, positions[:, 0])
            interp_y = CubicSpline(epochs_ref, positions[:, 1])
            interp_z = CubicSpline(epochs_ref, positions[:, 2])
            if np.sum(valid_clk) >= 2:
                from scipy.interpolate import interp1d

                interp_clock = interp1d(
                    epochs_ref[valid_clk],
                    positions[valid_clk, 3],
                    kind="linear",
                    bounds_error=False,
                    fill_value=np.nan,
                )
            else:
                interp_clock = None
        elif interp_type == "barycentric":
            interp_x = BarycentricInterpolator(epochs_ref, positions[:, 0])
            interp_y = BarycentricInterpolator(epochs_ref, positions[:, 1])
            interp_z = BarycentricInterpolator(epochs_ref, positions[:, 2])
            if np.sum(valid_clk) >= 2:
                from scipy.interpolate import interp1d

                interp_clock = interp1d(
                    epochs_ref[valid_clk],
                    positions[valid_clk, 3],
                    kind="linear",
                    bounds_error=False,
                    fill_value=np.nan,
                )
            else:
                interp_clock = None
        else:
            raise ValueError(f"Unsupported interpolation type: {interp_type}")

        clock = (
            interp_clock(epoch)
            if interp_clock is not None
            else np.full_like(np.asarray(epoch, dtype=float), np.nan)
        )

        # range of epochs
        min_epoch = np.min(epochs_ref)
        max_epoch = np.max(epochs_ref)

        if epoch_vec:
            if np.any(epoch < min_epoch) or np.any(epoch > max_epoch):
                propagate = True
        else:
            if epoch < min_epoch or epoch > max_epoch:
                propagate = True

        if propagate:
            if epoch_vec:
                # use the first epoch in the vector for propagation
                epoch0 = epoch[0]  # use the first epoch for propagation
            else:
                epoch0 = epoch

            idx = np.searchsorted(epochs_ref, epoch0)

            if idx == 0:
                r0_ecef = positions[0]
                t0 = epochs_ref[0]
            elif idx == len(epochs_ref):
                r0_ecef = positions[-1]
                t0 = epochs_ref[-1]
            else:
                r0_ecef = positions[idx - 1]
                t0 = epochs_ref[idx - 1]

            # compute velocity using
            if interp_type == "cubic":
                vx0 = interp_x.derivative()(t0)
                vy0 = interp_y.derivative()(t0)
                vz0 = interp_z.derivative()(t0)
            elif interp_type == "barycentric":
                vx0 = interp_x.derivative(t0)
                vy0 = interp_y.derivative(t0)
                vz0 = interp_z.derivative(t0)

            rv0_ecef = np.array(
                [r0_ecef[0], r0_ecef[1], r0_ecef[2], vx0, vy0, vz0]
            )  # ECEF position and velocity in km/s

            # convert to ECI frame
            t0_tdb = pnt.convert_time(t0, pnt.TAI, pnt.TDB)
            rv0_eci = pnt.convert_frame(t0_tdb, rv0_ecef, pnt.ECEF, pnt.ECI, rotate_only=False)

            # print("rv0_eci:", rv0_eci)
            # print(f"Using epoch {t0} for satellite {sat} at requested epoch {epoch}")

            # propagate the orbit to the requested epoch
            if not epoch_vec:
                if t0 == epoch:
                    rv_prop = rv0_eci
                else:
                    rv_prop = self.dyn_prop.propagate(rv0_eci, t0, epoch, np.zeros(6))
            else:
                rv_prop = self.dyn_prop.propagate(rv0_eci, t0, epoch)

            if out_frame != pnt.ECI:
                # convert to the requested frame
                epoch_tdb = pnt.convert_time(epoch, pnt.TAI, pnt.TDB)
                rv_prop = pnt.convert_frame(
                    epoch_tdb, rv_prop, pnt.ECI, out_frame, rotate_only=False
                )

        else:
            # use interpolation without propagation
            if interp_type == "cubic":
                vx = interp_x.derivative()(epoch)
                vy = interp_y.derivative()(epoch)
                vz = interp_z.derivative()(epoch)
            elif interp_type == "barycentric":
                vx = interp_x.derivative(epoch)
                vy = interp_y.derivative(epoch)
                vz = interp_z.derivative(epoch)

            rv_prop = np.vstack(
                (interp_x(epoch), interp_y(epoch), interp_z(epoch), vx, vy, vz)
            ).T  # ECEF position and velocity in km/s  T x 6

            if rv_prop.shape[0] == 1 and not epoch_vec:
                rv_prop = rv_prop[0]
            # if not epoch_vec:
            #     rv_prop = rv_prop[0]

            if out_frame != pnt.ECEF:
                # convert to the requested frame
                epoch_tdb = pnt.convert_time(epoch, pnt.TAI, pnt.TDB)
                rv_prop = pnt.convert_frame(
                    epoch_tdb, rv_prop, pnt.ECEF, out_frame, rotate_only=False
                )
                if rv_prop.shape[0] == 1 and not epoch_vec:
                    rv_prop = rv_prop[0]

        # add relativisic correction
        C_ms = 299792458  # speed of light in m/s

        if not epoch_vec:
            dot_rv = np.dot(rv_prop[:3], rv_prop[3:])  # dot product of position and velocity
        else:
            dot_rv = np.array([np.dot(rv[:3], rv[3:]) for rv in rv_prop])

        t_corr_rel = -2 / C_ms / C_ms * dot_rv  # in seconds

        clock = clock + t_corr_rel  # in seconds

        return rv_prop, clock

    def get_fault_clock_prns(self):
        """Return sorted list of SP3 satellite IDs that had at least one sentinel (999999.999999 µs) clock epoch."""
        return sorted(self.fault_clock_sats)

    def get_posvel(self, sys, prn, epoch, out_frame=pnt.ECI, propagate=False):
        """
        Get the position and velocity of a satellite at a specific epoch.
        """
        rv_prop, clock = self.get_posvelclock(
            sys, prn, epoch, out_frame=out_frame, propagate=propagate
        )
        return rv_prop

    def get_posvelclock_all(self, epochs_tai, out_frame=pnt.ECI, sats_list=None, propagate=False):
        """
        Get all satellite orbits for the given time span.

        returns:
        orbits: N x T x 6 array of satellite positions and velocities
        clocks: N x T array of clock corrections
        """
        if sats_list is None:
            sats_list = [sat for sat in self.sats]  # exclude GLONASS satellites for now

        nsats = len(sats_list)
        lent = epochs_tai.shape[0]
        orbits = np.zeros((nsats, lent, 6))  # N x T x 6 array of satellite positions and velocities
        clocks = np.zeros((nsats, lent))  # N x T array of clock corrections

        for i, sat in enumerate(sats_list):
            orbits[i], clocks[i] = self.get_posvelclock(
                sat[0], int(sat[1:]), epochs_tai, out_frame=out_frame, propagate=propagate
            )

        return orbits, clocks, sats_list


class BRDCLoader:
    """
    Class to load and parse Brodcast ephemeris files.
    """

    def __init__(self, filenames=[], target_dt=None, sim_t=0, dt_timesys=pnt.UTC):
        """
        Initialize the BRDCLoader with the filename of the BRDC file.
        """
        self.datapath = os.path.join(pnt.get_output_dir("gnss_files"), "brdc")
        if not os.path.exists(self.datapath):
            os.makedirs(self.datapath)

        self.first_load = True  # flag to indicate if this is the first load

        is_list_filenames = isinstance(filenames, list)

        if not is_list_filenames:
            if filenames is None:
                filenames = []
            if isinstance(filenames, str):
                filenames = [filenames]

        if sim_t == 0:
            target_dts = [target_dt] if target_dt is not None else []
        else:
            target_dt_end = (
                target_dt + timedelta(seconds=sim_t) + timedelta(days=1)
            )  # margin of 1 day
            target_dt_start = target_dt - timedelta(days=1)  # margin of 1 day

            target_dt = target_dt_start
            target_dts = []
            while target_dt <= target_dt_end:
                target_dts.append(target_dt)
                target_dt += timedelta(days=1)

        if len(filenames) == 0 and len(target_dts) == 0:
            # error if no filename or date is provided
            raise ValueError("Either filename or date must be provided.")
        elif len(filenames) == 0 and len(target_dts) > 0:
            # construct filename from date
            self.filenames = [self.load_brdc(dt, dt_timesys) for dt in target_dts]
        elif len(filenames) > 0:
            self.filenames = filenames
        else:
            raise ValueError("Invalid Input")

        # parse the BRDC file
        self.iono_dict = {}
        self.nav_dict = {}
        self.sats = []  # use a set to avoid duplicates

        for filename in self.filenames:
            self.parse_brdc(filename)

    def load_brdc(self, target_dt, dt_timesys=pnt.UTC):
        """
        Load the BRDC file from the CDDIS archive. (RINEX V3 format)
        Reference: https://www.earthdata.nasa.gov/data/space-geodesy-techniques/gnss/broadcast-ephemeris-data-product
        """

        _, _, _, target_dt = datetime_to_gpsweeks(target_dt, dt_timesys)  # convert to GPS time

        # Construct the filename based on the date
        year = target_dt.year
        month = target_dt.month
        doy = target_dt.timetuple().tm_yday

        YYYY = str(year)
        YY = YYYY[2:4]  # last two digits of the year
        DDD = f"{doy:03d}"  # day of the year, zero-padded to 3 digits

        filename = f"BRDC00IGS_R_{YYYY}{DDD}0000_01D_MN.rnx.gz"
        filename1 = f"{YYYY}/{DDD}/{YY}p/{filename}"
        filename2 = f"{YYYY}/brdc/{filename}"

        if os.path.exists(os.path.join(self.datapath, filename[:-3])):
            print(
                f"BRDC file for GPS TIME: {target_dt} already exists at {os.path.join(self.datapath, filename[:-3])}. Skipping download."
            )
            return os.path.join(self.datapath, filename[:-3])

        # Download the file from CDDIS
        url = f"https://cddis.nasa.gov/archive/gnss/data/daily/{filename1}"
        response = requests.get(url)

        if response.status_code == 200:
            filesave = os.path.join(self.datapath, filename)
            with open(filesave, "wb") as f:
                f.write(response.content)
            # unzpip the file if it is gzipped
            if filesave.endswith(".gz"):
                with gzip.open(filesave, "rb") as f_in:
                    with open(filesave[:-3], "wb") as f_out:
                        f_out.write(f_in.read())
                os.remove(filesave)
                filesave = filesave[:-3]  # Remove .gz extension
            print(f"BRDC file downloaded from {url} and saved as {filesave}")
            return filesave
        else:
            # try second URL
            url2 = f"https://cddis.nasa.gov/archive/gnss/data/daily/{filename2}"
            response = requests.get(url2)
            filesave = os.path.join(self.datapath, filename)
            if response.status_code == 200:
                with open(filesave, "wb") as f:
                    f.write(response.content)

                # unzpip the file if it is gzipped
                if filesave.endswith(".gz"):
                    with gzip.open(filesave, "rb") as f_in:
                        with open(filesave[:-3], "wb") as f_out:
                            f_out.write(f_in.read())
                    os.remove(filesave)
                    filesave = filesave[:-3]  # Remove .gz extension
                print(f"BRDC file downloaded from {url2} and saved as {filesave}")
                return filesave
            else:
                raise FileNotFoundError(
                    f"Brodcast file for {target_dt}(url:{url} or {url2}) not found in CDDIS archive."
                )

    def parse_brdc(self, filename):

        print(f"Parsing BRDC file: {filename}")

        # Named orbit labels by system
        # http://acc.igs.org/misc/rinex304.pdf
        ORBITS = {
            # GPS (A6)
            "G": [
                "epoch_dt",
                "epoch_tai",
                "af0",
                "af1",
                "af2",
                "IODE",
                "Crs",
                "Delta_n",
                "M0",
                "Cuc",
                "ecc",
                "Cus",
                "sqrtA",
                "Toe",
                "Cic",
                "Omega0",
                "Cis",
                "i0",
                "Crc",
                "omega",
                "Omega_dot",
                "IDOT",
                "Codes_L2",
                "Week",
                "L2P_flag",
                "SV_accuracy",
                "SV_health",
                "TGD",
                "IODC",
                "T_trans",
                "Fit",
            ],
            # Galileo (A8)
            "E": [
                "epoch_dt",
                "epoch_tai",
                "af0",
                "af1",
                "af2",
                "IODnav",
                "Crs",
                "Delta_n",
                "M0",
                "Cuc",
                "ecc",
                "Cus",
                "sqrtA",
                "Toe",
                "Cic",
                "Omega0",
                "Cis",
                "i0",
                "Crc",
                "omega",
                "Omega_dot",
                "IDOT",
                "DataSrc",
                "Week",
                "SV_accuracy",
                "SV_health",
                "BGD_E5aE1",
                "BGD_E5bE1",
                "T_trans",
            ],
            # GLONASS (A10)
            "R": [
                "epoch_dt",
                "epoch_tai",
                "tau_n",
                "gamma_n",
                "tk",
                "X",
                "VX",
                "AX",
                "health",
                "Y",
                "VY",
                "AY",
                "freq",
                "Z",
                "VZ",
                "AZ",
                "age",
            ],
            # BeiDou (A14)
            "C": [
                "epoch_dt",
                "epoch_tai",
                "af0",
                "af1",
                "af2",
                "AODE",
                "Crs",
                "Delta_n",
                "M0",
                "Cuc",
                "ecc",
                "Cus",
                "sqrtA",
                "Toe",
                "Cic",
                "Omega0",
                "Cis",
                "i0",
                "Crc",
                "omega",
                "Omega_dot",
                "IDOT",
                "Spare1",
                "Week",
                "Spare2",
                "SV_accuracy",
                "SV_health",
                "TGD_B1B3",
                "TGD_B2B3",
                "T_trans",
                "AODC",
            ],
            # QZSS (A12)
            "J": [
                "epoch_dt",
                "epoch_tai",
                "af0",
                "af1",
                "af2",
                "IODE",
                "Crs",
                "Delta_n",
                "M0",
                "Cuc",
                "ecc",
                "Cus",
                "sqrtA",
                "Toe",
                "Cic",
                "Omega0",
                "Cis",
                "i0",
                "Crc",
                "omega",
                "Omega_dot",
                "IDOT",
                "Codes_L2",
                "Week",
                "L2P_flag",
                "SV_accuracy",
                "SV_health",
                "TGD",
                "IODC",
                "T_trans",
                "Fit",
            ],
        }

        nav = {"G": {}, "E": {}, "C": {}, "J": {}, "R": {}}

        iono = {}
        times = {}

        with open(filename) as f:
            lines = f.readlines()

        # parse the ionospheric correction coefficients and leap seconds ---------------
        # print("Parsing ionospheric correction coefficients and leap seconds")
        i = 0
        while i < len(lines) and "END OF HEADER" not in lines[i]:
            line = lines[i]
            # parse the ionospheric coefficients
            if "IONOSPHERIC" in line:
                # split line by whitespace and extract coefficients
                parts = line.split()

                iono[parts[0]] = []

                for j in range(1, len(parts)):
                    # if convertable to float, add to iono
                    try:
                        iono[parts[0]].append(float(parts[j]))
                    except ValueError:
                        break  # stop if not a float

            if "TIME SYSTEM CORR" in line:
                FLOAT_RE = re.compile(
                    r"""
                        [+-]?                 # optional sign
                        (?:\d+\.\d*|\.\d+|\d+) # mantissa (supports 1., .1, 1.23, 123)
                        (?:[Ee][+-]?\d+)?     # optional exponent
                    """,
                    re.VERBOSE,
                )

                parts = line.split()
                key = parts[0]
                times[key] = []

                # Only parse the section before the label (safer than relying on split positions)
                left = line.split("TIME SYSTEM CORR", 1)[0]
                left = left.replace("D", "E").replace("d", "e")

                nums = [float(x) for x in FLOAT_RE.findall(left)]

                # In your example this yields: [0.0, -9.769962617e-15, 604385.0, 1000.0]
                # If you only want the *corrections* (usually first two):
                times[key].extend(nums)

                # # split line by whitespace and extract time system corrections
                # parts = line.split()
                # times[parts[0]] = []

                # for j in range(1, len(parts) - 1):
                #     # if length of the part is larger than 12, two minus signs are connected
                #     if len(parts[j]) > 20 and "-" in parts[j]:
                #         subparts = parts[j].split("-")
                #         n_nums = 0
                #         for k in range(len(subparts)):
                #             if ("D" in subparts[k]) or ("d" in subparts[k]):
                #                 # replace D with E for scientific notation
                #                 subparts[k] = subparts[k].replace("D", "E").replace("d", "e")
                #             if ("E" in subparts[k]) or ("e" in subparts[k]):
                #                 if k == 0:  # the first number is positive
                #                     num_parts = (
                #                         subparts[k] + "-" + subparts[k + 1]
                #                     )  # combine with the next part
                #                 else:  # the number is negative
                #                     num_parts = (
                #                         "-" + subparts[k] + "-" + subparts[k + 1]
                #                     )  # combine with the next part
                #                 # convert scientific notation to float
                #                 times[parts[0]].append(float(num_parts))

                #     # if convertable to float, add to times
                #     else:
                #         try:
                #             times[parts[0]].append(float(parts[j]))
                #         except ValueError:
                #             break

            # if "LEAP SECONDS" in line:
            #     leap = int(line.split()[0])

            i += 1

        # summarize results
        i += 1
        # print("Leap seconds:", leap)
        # print("Ionospheric coefficients:", iono)

        # parse the ephemeris data ------------------------------------------------------
        # print("Parsing ephemeris data")
        recorded_ymd = False

        while i < len(lines):
            line = lines[i]
            init_str = line[0]
            if init_str not in ("G", "E", "C", "J", "R"):
                i += 1
                continue

            # split by space
            parts = line.split()
            sys = parts[0][0]  # system identifier
            prn = int(parts[0][1:])  # PRN number
            year = int(parts[1])
            month = int(parts[2])
            day = int(parts[3])
            hh = int(parts[4])
            mm = int(parts[5])

            if not recorded_ymd:
                sim_year = year
                sim_month = month
                sim_day = day
                sim_hour = hh
                sim_min = mm
                recorded_ymd = True

            # second can be connected to the next text if it is not separated by space
            #  G02 2020 10 01 00 00 00-5.236091092229E-04-5.115907697473E-12 0.000000000000E+00
            ss = int(line[20:22])
            af0 = float(line[23:42])  # af0
            af1 = float(line[42:61])  # af1
            af2 = float(line[61:80])  # af2

            epoch_dt = datetime(year, month, day, hh, mm, ss)
            epoch_tai = datetime_to_tai(epoch_dt, dt_timesys=pnt.Time.GPS)
            # print("Parsing data for system:", sys, "PRN:", prn, "Epoch:", epoch)

            if sys in ("G", "E", "C", "J"):
                vals = []
                for k in range(7):
                    i += 1
                    line = lines[i].strip()
                    start = 0
                    for kk in range(4):
                        if kk == 0 and line[0] == "-":
                            lennum = 19  # consider the sign
                        elif kk == 0:
                            lennum = 18
                        else:
                            lennum = 19

                        val = line[start : start + lennum].strip()
                        if val:
                            vals.append(float(val))
                        start = start + lennum

                fields = ORBITS[sys][5:]
                rec = {
                    "epoch_dt": epoch_dt,
                    "epoch_tai": epoch_tai,
                    "af0": af0,
                    "af1": af1,
                    "af2": af2,
                }
                for field, val in zip(fields, vals):
                    rec[field] = val

                if prn not in nav[sys]:
                    nav[sys][prn] = {}
                    for key in ORBITS[sys]:
                        nav[sys][prn][key] = []

                for key, value in rec.items():
                    nav[sys][prn][key].append(value)

                i += 1

            elif sys == "R":
                vals = []
                for k in range(3):
                    i += 1
                    line = lines[i].strip()
                    start = 0
                    for kk in range(4):
                        if kk == 0 and line[0] == "-":
                            lennum = 19  # consider the sign
                        elif kk == 0:
                            lennum = 18
                        else:
                            lennum = 19

                        # if len(line) > start + lennum:
                        val = line[start : start + lennum].strip()
                        if val:
                            vals.append(float(val))
                        start = start + lennum

                fields = ORBITS[sys][5:]
                rec = {
                    "epoch_dt": epoch_dt,
                    "epoch_tai": epoch_tai,
                    "tau_n": af0,
                    "gamma_n": af1,
                    "tk": af2,
                }
                for field, val in zip(fields[3:], vals):
                    rec[field] = val

                if prn not in nav[sys]:
                    nav[sys][prn] = {}
                    for key in ORBITS[sys]:
                        nav[sys][prn][key] = []

                for key, value in rec.items():
                    nav[sys][prn][key].append(value)
                i += 1
            else:
                i += 1

        # add year, month, day to the ionospheric coefficients
        iono["year"] = sim_year
        iono["month"] = sim_month
        iono["day"] = sim_day
        times["year"] = sim_year
        times["month"] = sim_month
        times["day"] = sim_day
        times["hr"] = sim_hour
        times["min"] = sim_min

        # convert the lists to numpy arrays for each satellite
        for sys in nav:
            for prn in nav[sys]:
                for key in ORBITS[sys]:
                    if nav[sys][prn][key]:
                        nav[sys][prn][key] = np.array(nav[sys][prn][key])
                    else:
                        nav[sys][prn][key] = np.array([])

        if self.first_load:
            self.iono_dicts = {}
            for key in iono:
                self.iono_dict[key] = [iono[key]]
            self.time_dicts = {}
            for key in times:
                self.time_dicts[key] = [times[key]]

            self.nav_dict = nav
            self.first_load = False
        else:
            # merge the dictionaries
            for key in iono:
                if key not in self.iono_dict:
                    self.iono_dict[key] = iono[key]
                else:
                    self.iono_dict[key].append(iono[key])

            for key in times:
                if key not in self.time_dicts:
                    self.time_dicts[key] = times[key]
                else:
                    self.time_dicts[key].append(times[key])

            for sys in nav:
                if sys not in self.nav_dict:
                    self.nav_dict[sys] = nav[sys]
                else:
                    for prn in nav[sys]:
                        if prn not in self.nav_dict[sys]:
                            self.nav_dict[sys][prn] = nav[sys][prn]
                        else:
                            epochs_nav = self.nav_dict[sys][prn]["epoch_tai"]
                            epochs_new = nav[sys][prn]["epoch_tai"]
                            # check if the epochs are already present
                            idx_new = [i for i, e in enumerate(epochs_new) if e not in epochs_nav]

                            for key in ORBITS[sys]:
                                if key not in self.nav_dict[sys][prn]:
                                    self.nav_dict[sys][prn][key] = nav[sys][prn][key]
                                else:
                                    idx_new = [
                                        i for i in idx_new if i < nav[sys][prn][key].shape[0]
                                    ]
                                    self.nav_dict[sys][prn][key] = np.concatenate(
                                        (self.nav_dict[sys][prn][key], nav[sys][prn][key][idx_new])
                                    )
        # constrct a list of satellites
        self.sats = []
        for sys in self.nav_dict:
            for prn in self.nav_dict[sys]:
                self.sats.append(f"{sys}{prn:02d}")

    def get_posvelclock(self, sys, prn, epoch_tai, out_frame=pnt.ECI):
        """
        Get the navigation data for a specific satellite at a specific epoch.

        Parameters:
        sys (str): Satellite system identifier (e.g., 'G' for GPS, 'E' for Galileo).
        prn (int): PRN number of the satellite.
        epoch (datetime): Epoch in UTC datetime format.

        Returns:
        dict: Navigation data for the satellite at the specified epoch.
        """

        if isinstance(epoch_tai, (list, np.ndarray)):
            epoch_vec = True
            if isinstance(epoch_tai, list):
                epoch_tai = np.array(epoch_tai)
        else:
            epoch_vec = False

        # constants used
        GM_EARTH = 3.986005e14  # m^3/s^2
        OMEGA_DOT_EARTH = 7.2921151467e-5  # rad/s, Earth's rotation rate

        if sys not in self.nav_dict or prn not in self.nav_dict[sys]:
            raise ValueError(f"Satellite {sys}{prn} not found in navigation data.")

        nav_data = self.nav_dict[sys][prn]

        # compute gps week and seconds into week
        gps_week, sec_week = tai_to_gps_weeks(epoch_tai)
        # print(f"Epoch TAI: {epoch_tai}, GPS Week: {gps_week}, Seconds into week: {sec_week}")

        # compute the year, month, day
        if epoch_vec:
            dt_epochs = [tai_to_datetime(et, dt_timesys=pnt.Time.GPS) for et in epoch_tai]
        else:
            dt_epochs = tai_to_datetime(epoch_tai, dt_timesys=pnt.Time.GPS)

        # Find the closest epoch
        if not epoch_vec:
            idx = np.argmin(np.abs(nav_data["epoch_tai"] - epoch_tai))
        else:
            idx = []
            for et in epoch_tai:
                idx.append(np.argmin(np.abs(nav_data["epoch_tai"] - et)))

        # epoch_dt = tai_to_datetime(epoch_tai, dt_timesys=pnt.Time.GPS)
        # print(f"Using GPS epoch {nav_data['epoch_dt'][idx]} for satellite {sys}{prn} at requested GPS epoch {epoch_dt}")

        # Apply keplarian transformation to get position and velocity
        # Extract the necessary parameters
        if sys in {"G", "E", "C", "J"}:
            # Extract parameters for GPS, Galileo, BeiDou, and QZSS
            af0 = nav_data["af0"][idx]
            af1 = nav_data["af1"][idx]
            af2 = nav_data["af2"][idx]

            sqrt_sma = nav_data["sqrtA"][idx]
            ecc = nav_data["ecc"][idx]
            i0 = nav_data["i0"][idx]
            Omega0 = nav_data["Omega0"][idx]
            omega = nav_data["omega"][idx]
            M0 = nav_data["M0"][idx]

            delta_n = nav_data["Delta_n"][idx]
            idot = nav_data["IDOT"][idx]
            Omega_dot = nav_data["Omega_dot"][idx]

            crs = nav_data["Crs"][idx]
            crc = nav_data["Crc"][idx]
            cus = nav_data["Cus"][idx]
            cuc = nav_data["Cuc"][idx]
            cis = nav_data["Cis"][idx]
            cic = nav_data["Cic"][idx]

            # Compute the time offset from the reference epoch
            t_oe = nav_data["Toe"][idx]  # sec of gps week
            t_k = sec_week - t_oe
            # print(f"t_k: {t_k} seconds")

            if epoch_vec:
                t_k[t_k > 302400] -= 604800  # if t_k is greater than 4 days, subtract 1 week
                t_k[t_k < -302400] += 604800  # if t_k is less than -4 days, add 1 week
            else:
                if t_k > 302400:  # if t_k is greater than 4 days, subtract 604800 (1 week)
                    t_k -= 604800
                elif t_k < -302400:  # if t_k is less than -4 days, add 604800 (1 week)
                    t_k += 604800

            # Anomaly calculations
            n0 = np.sqrt(GM_EARTH) * np.ones_like(sqrt_sma) / (sqrt_sma**3)  # Mean motion
            n = n0 + delta_n  # Corrected mean motion
            M = M0 + n * t_k  # Mean anomaly

            # Solve Kepler's equation for eccentric anomaly
            E = M
            for _ in range(3):  # Iterate to solve for E
                E = E + (M - (E - ecc * np.sin(E))) / (1 - ecc * np.cos(E))
            # True anomaly
            nu = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E / 2))

            # Time correction
            t_corr_poly = af0 + af1 * t_k + af2 * t_k**2

            # Relativistic correction (not needed to match with use)
            F = -4.442807633e-10  # relativistic correction factor (s/sqrt(m))
            t_corr_rel = F * ecc * sqrt_sma * np.sin(E)

            if sys == "E":
                # GALILEO specific time corrections
                if isinstance(idx, (list, np.ndarray)):
                    time_dicts_idx = []
                    num_idx = len(idx)
                    num_time_dicts = len(self.time_dicts["year"])
                    for ii in range(num_idx):
                        time_search = datetime(
                            dt_epochs[ii].year, dt_epochs[ii].month, dt_epochs[ii].day, 0, 0, 0
                        )
                        for kk in range(num_time_dicts):
                            time_query = datetime(
                                self.time_dicts["year"][kk],
                                self.time_dicts["month"][kk],
                                self.time_dicts["day"][kk],
                                self.time_dicts["hr"][kk],
                                self.time_dicts["min"][kk],
                                0,
                            )
                            # find the index of the time dicts where the year, month, day match
                            time_diff = abs((time_query - time_search).total_seconds())
                            if time_diff < 1 * 3600:  # within 1 hour
                                time_dicts_idx.append(kk)
                                break
                    coeffs = np.zeros(
                        (len(time_dicts_idx), len(self.time_dicts["GAGP"][0]))
                    )  # coefficients for each time dict
                    for kk, tidx in enumerate(time_dicts_idx):
                        coeffs[kk, :] = np.array(self.time_dicts["GAGP"])[tidx]
                    t_corr_sys = coeffs[:, 0] + coeffs[:, 1] * (
                        (sec_week - coeffs[:, 2]) + 604800 * (gps_week - coeffs[:, 3])
                    )
                else:
                    years = dt_epochs.year
                    months = dt_epochs.month
                    days = dt_epochs.day
                    time_search = datetime(years, months, days, 0, 0, 0)
                    # print(f"Years: {years}, Months: {months}, Days: {days}")
                    year_td = self.time_dicts["year"]
                    month_td = self.time_dicts["month"]
                    day_td = self.time_dicts["day"]
                    hh_td = self.time_dicts["hr"]
                    min_td = self.time_dicts["min"]
                    # print(f"Time dicts length: {len(year_td)}"  )
                    for ii in range(len(year_td)):
                        time_query = datetime(
                            year_td[ii], month_td[ii], day_td[ii], hh_td[ii], min_td[ii], 0
                        )
                        time_diff = abs((time_query - time_search).total_seconds())
                        # print(f"Comparing to time dict {ii}: {time_query}, diff: {time_diff} seconds")
                        if time_diff < 1 * 3600:  # within 1 hour
                            time_dicts_idx = ii
                            break
                    # print(f"Using time dict index: {time_dicts_idx}")
                    coeffs = np.array(self.time_dicts["GAGP"][time_dicts_idx])
                    # print(f"Using coefficients: {coeffs}")
                    t_corr_sys = coeffs[0] + coeffs[1] * (
                        (sec_week - coeffs[2]) + 604800 * (gps_week - coeffs[3])
                    )

            else:
                # For GPS, BeiDou, and QZSS, we use the same correction as in the original code
                t_corr_sys = 0

            # t_corr_rel = 0
            t_clk_corr = t_corr_poly + t_corr_rel + t_corr_sys  # in seconds

            # Calculate the satellite position in the orbital plane
            A = sqrt_sma**2
            phi = nu + omega  # Argument of latitude

            delta_u = cus * np.sin(2 * phi) + cuc * np.cos(2 * phi)
            delta_r = crs * np.sin(2 * phi) + crc * np.cos(2 * phi)
            delta_i = cis * np.sin(2 * phi) + cic * np.cos(2 * phi)

            u_k = phi + delta_u
            r_k = A * (1 - ecc * np.cos(E)) + delta_r
            i_k = i0 + delta_i + idot * t_k
            Omega_k = Omega0 + (Omega_dot - OMEGA_DOT_EARTH) * t_k - OMEGA_DOT_EARTH * t_oe

            x_k_hat = r_k * np.cos(u_k)
            y_k_hat = r_k * np.sin(u_k)

            # positions
            x_k = x_k_hat * np.cos(Omega_k) - y_k_hat * np.cos(i_k) * np.sin(Omega_k)
            y_k = x_k_hat * np.sin(Omega_k) + y_k_hat * np.cos(i_k) * np.cos(Omega_k)
            z_k = y_k_hat * np.sin(i_k)

            # velocities
            Ekdot = n / (1 - ecc * np.cos(E))  # derivative of E with respect to time
            nudot = (
                Ekdot * np.sqrt(1 - ecc**2) / (1 - ecc * np.cos(E))
            )  # derivative of nu with respect to time
            didot_dt = idot + 2 * nudot * (
                cis * np.cos(2 * phi) - cic * np.sin(2 * phi)
            )  # derivative of inclination with respect to time
            udot = nudot + 2 * nudot * (
                cus * np.cos(2 * phi) - cuc * np.sin(2 * phi)
            )  # derivative of argument of latitude with respect to time
            rdot = ecc * A * Ekdot * np.sin(E) + 2 * nudot * (
                crs * np.cos(2 * phi) - crc * np.sin(2 * phi)
            )  # derivative of radius with respect to time
            Omega_k_dot = (
                Omega_dot - OMEGA_DOT_EARTH
            )  # derivative of right ascension with respect to time

            xdot_hat = rdot * np.cos(u_k) - r_k * udot * np.sin(u_k)
            ydot_hat = rdot * np.sin(u_k) + r_k * udot * np.cos(u_k)

            xdot_k = (
                -x_k_hat * Omega_k_dot * np.sin(Omega_k)
                + xdot_hat * np.cos(Omega_k)
                - ydot_hat * np.sin(Omega_k) * np.cos(i_k)
                - y_k_hat
                * (
                    Omega_k_dot * np.cos(Omega_k) * np.cos(i_k)
                    - didot_dt * np.sin(Omega_k) * np.sin(i_k)
                )
            )
            ydot_k = (
                x_k_hat * Omega_k_dot * np.cos(Omega_k)
                + xdot_hat * np.sin(Omega_k)
                + ydot_hat * np.cos(Omega_k) * np.cos(i_k)
                - y_k_hat
                * (
                    Omega_k_dot * np.sin(Omega_k) * np.cos(i_k)
                    + didot_dt * np.cos(Omega_k) * np.sin(i_k)
                )
            )
            zdot_k = y_k_hat * didot_dt * np.cos(i_k) + ydot_hat * np.sin(i_k)

            # concatenate position and velocity
            rv_ecef = np.vstack([x_k, y_k, z_k, xdot_k, ydot_k, zdot_k]).T  # convert to m/s (T x 6)

            if not epoch_vec:
                rv_ecef = rv_ecef[0]

            if out_frame != pnt.ECEF:
                # convert to the requested frame
                epoch_tdb = pnt.convert_time(epoch_tai, pnt.TAI, pnt.TDB)
                rv_out = pnt.convert_frame(
                    epoch_tdb, rv_ecef, pnt.ECEF, out_frame, rotate_only=False
                )
                if not epoch_vec:
                    rv_out = rv_out[0]
            else:
                rv_out = rv_ecef

            return rv_out, t_clk_corr

        else:
            raise ValueError(f"Satellite system {sys} not supported.")

    def get_posvel(self, sys, prn, epoch, out_frame=pnt.ECI):
        """
        Get the position and velocity of a satellite at a specific epoch.
        """
        rv_prop, clock = self.get_posvelclock(sys, prn, epoch, out_frame=out_frame)
        return rv_prop

    def get_posvelclock_all(self, epochs_tai, out_frame=pnt.ECI, sats_list=None):
        """
        Get all satellite orbits for the given time span.

        returns:
        orbits: N x T x 6 array of satellite positions and velocities
        clocks: N x T array of clock corrections
        """
        if sats_list is None:
            sats_list = [
                sat for sat in self.sats if sat[0] != "R"
            ]  # exclude GLONASS satellites for now

        nsats = len(sats_list)
        lent = epochs_tai.shape[0]
        orbits = np.zeros((nsats, lent, 6))  # N x T x 6 array of satellite positions and velocities
        clocks = np.zeros((nsats, lent))  # N x T array of clock corrections

        for i, sat in enumerate(sats_list):
            orbits[i], clocks[i] = self.get_posvelclock(
                sat[0], int(sat[1:]), epochs_tai, out_frame=out_frame
            )

        return orbits, clocks, sats_list
