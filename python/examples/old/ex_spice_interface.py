"""Call the public SPICE-backed helpers exposed through ``pylupnt``.

The direct ``pylupnt.spice`` submodule is not exported by the current package,
so this example uses the supported top-level time, frame, and ephemeris APIs.
"""

import pylupnt as pnt


def main() -> None:
    # Use the numeric Gregorian overload; it round-trips cleanly in current builds.
    t_tai = pnt.gregorian_to_time(2000, 1, 1, 12, 0, 0.0)
    t_tdb = pnt.convert_time(t_tai, pnt.Time.TAI, pnt.Time.TDB)
    utc = pnt.time_to_gregorian_string(t_tai, 3)

    # Body states and frame transforms use TDB coordinate time.
    earth_wrt_sun = pnt.get_body_pos_vel(t_tdb, pnt.SUN, pnt.EARTH, pnt.Frame.ICRF)
    rot_gcrf_to_itrf, offset = pnt.get_frame_rotation_translation(
        t_tdb, pnt.Frame.GCRF, pnt.Frame.ITRF
    )

    print("TDB at 2000/01/01 12:00:00  :", t_tdb)
    print("TAI at 2000/01/01 12:00:00  :", t_tai)
    print("UTC at above TAI            :", utc)
    print("Earth relative to Sun       :")
    print(earth_wrt_sun)
    print("Rotation from GCRF to ITRF  :")
    print(rot_gcrf_to_itrf)
    print("Translation from GCRF to ITRF:")
    print(offset)


if __name__ == "__main__":
    main()
