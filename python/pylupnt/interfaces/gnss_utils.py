from datetime import datetime
import pylupnt as pnt


def time_to_mjd(t):
    SECS_DAY = 86400.0  # seconds in a day
    MJD_J2000_TT = 51544.5
    # [days]
    return t / SECS_DAY + MJD_J2000_TT


def datetime_to_gpsweeks(dt: datetime, timesys=pnt.UTC):
    """
    Convert UTC datetime to GPS week number and seconds into week.

    :returns: ``(week, day_of_week, sec_of_week, dt_gps)`` where ``week`` is
        the GPS week number, ``day_of_week`` is 0 for Sunday through 6 for
        Saturday, ``sec_of_week`` is the seconds into the week, and ``dt_gps``
        is the corresponding datetime in the GPS time system.
    """

    # GPS epoch: Sunday, Jan 6 1980
    gps_epoch_dt = datetime(1980, 1, 6, 0, 0, 0)
    gps_epoch_tai = datetime_to_tai(gps_epoch_dt, dt_timesys=pnt.Time.GPS)

    curr_tai = datetime_to_tai(dt, dt_timesys=timesys)

    # Align dt to GPS time by adding leap seconds
    delta = curr_tai - gps_epoch_tai
    delta_days = delta // 86400

    # Integer division for weeks and seconds
    week = int(delta_days // 7)
    day_of_week = int(delta_days % 7)  # 0 = Sunday, 1 = Monday, ..., 6 = Saturday
    sec_of_week = delta - (week * 7 * 86400)

    dt_gps = tai_to_datetime(curr_tai, dt_timesys=pnt.Time.GPS)

    return week, day_of_week, sec_of_week, dt_gps


def gps_weeks_to_tai(week, sec_into_week):
    """
    Convert GPS week number and seconds into week to TAI seconds.
    Returns:
        tai_seconds (float): TAI seconds
    """
    gps_epoch_dt = datetime(1980, 1, 6, 0, 0, 0)
    gps_epoch_tai = datetime_to_tai(gps_epoch_dt, dt_timesys=pnt.Time.GPS)

    tai_seconds = gps_epoch_tai + week * 7 * 86400 + sec_into_week

    return tai_seconds


def gps_weeks_to_datetime(week, sec_into_week, dt_timesys=pnt.UTC):
    """
    Convert GPS week number and seconds into week to UTC datetime.
    Returns:
        dt (datetime): Corresponding datetime in specified time system
    """
    tai_seconds = gps_weeks_to_tai(week, sec_into_week)
    dt = tai_to_datetime(tai_seconds, dt_timesys=dt_timesys)
    return dt


def tai_to_gps_weeks(tai_seconds):
    """
    Convert TAI seconds to GPS week number and seconds into week.

    Returns:
        week (int): GPS week number
        sec_into_week (float): Seconds into the week
    """
    gps_epoch_dt = datetime(1980, 1, 6, 0, 0, 0)
    gps_epoch_tai = datetime_to_tai(gps_epoch_dt, dt_timesys=pnt.Time.GPS)

    delta = tai_seconds - gps_epoch_tai

    # Calculate the number of days since the GPS epoch
    delta_days = delta // 86400

    # Integer division for weeks and seconds
    week = delta_days // 7
    sec_into_week = delta - (week * 7 * 86400)  # seconds into the week

    return week, sec_into_week


def datetime_to_tai(dt: datetime, dt_timesys=pnt.UTC):
    """
    Convert UTC datetime to TAI seconds.
    """
    year = dt.year
    month = dt.month
    day = dt.day
    hour = dt.hour
    minute = dt.minute
    second = dt.second + dt.microsecond / 1e6  # Convert microseconds

    epoch = pnt.gregorian_to_time(year, month, day, hour, minute, second)

    if dt_timesys != pnt.TAI:
        # If the input time system is not TAI, convert to TAI
        epoch = pnt.convert_time(epoch, dt_timesys, pnt.TAI)

    return epoch


def tai_to_datetime(tai_seconds, dt_timesys=pnt.UTC):
    """
    Convert TAI seconds to UTC datetime.
    """
    # Convert TAI seconds to Gregorian time
    if dt_timesys != pnt.TAI:
        epoch = pnt.convert_time(tai_seconds, pnt.TAI, dt_timesys)

    # convert to mjd
    mjd = time_to_mjd(epoch)

    # Convert MJD to datetime
    dt_tuple = pnt.mjd_to_gregorian(mjd)

    year, month, day, hour, minute, second = dt_tuple

    # Create datetime object
    dt = datetime(year, month, day, hour, minute, int(second), int((second - int(second)) * 1e6))

    return dt
