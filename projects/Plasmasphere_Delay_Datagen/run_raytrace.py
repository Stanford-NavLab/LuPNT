import numpy as np
import pylupnt as pnt
import tecsimpy as tec
import multiprocessing as mp
import os


def run_single_job(
    job_id,
    job_nums,
    epoch_str,
    sat_idx,
    gnss_const,
    freq_family,
    Rz12,
    kp,
    correction,
    straight_ray,
    data_path,
    save_path,
):
    """
    Wrapper that calls tecsim.run_raytrace_batch for a single job_id.
    """
    print(
        f"[job {job_id}] Starting: job_id={job_id}, job_nums={job_nums}, "
        f"epoch={epoch_str}, sat_idx={sat_idx}, gnss_const={gnss_const}, "
        f"freq_family={freq_family}, Rz12={Rz12}, kp={kp}, "
        f"correction={correction}, straight_ray={straight_ray}",
        flush=True,
    )

    try:
        tec.run_raytrace_batch(
            job_id,
            job_nums,
            epoch_str,
            sat_idx,
            gnss_const,
            freq_family,
            Rz12,
            kp,
            correction,
            straight_ray,
            data_path,
            save_path,
        )
    except Exception as e:
        print(f"[job {job_id}] FAILED with exception: {e}", flush=True)
        # re-raise so the parent can notice
        raise

    print(f"[job {job_id}] Finished successfully", flush=True)
    return job_id


if __name__ == "__main__":
    # Inputs
    epoch_str = "2025_03_01_12_00_00"  # Epoch in UTC
    job_nums = 10  # Number of parallel jobs
    sat_nums = [0]  # Satellite numbers to process
    gnss_consts = [
        "GPS",
        "GALILEO",
        "QZSS",
    ]  # , "GALILEO", "QZSS"]  # GNSS constellation ("GPS", "Galileo", etc.)
    freq_families = [1]  # Frequency family ("L1"=1, "L2"=2, "L5"=5)
    correction = False  # Whether to apply correction (Always set to False)
    parallel = True  # Whether to run in parallel

    # Setting this to False makes you the simulation run for only min_alt < 2000 km
    straight_rays = [
        True
    ]  # Whether to assume straight ray path  <========================= CHANGE HERE
    Rz12s = [150.0]  # Rz12 index for ionospheric model
    kp = 3.0  # Kp index for geomagnetic activity

    n_orbit = 6
    dt = 1
    dtrt = 120

    # On macOS / recent Python, use 'spawn' to be safe with C++/Fortran libs
    try:
        mp.set_start_method("spawn")
    except RuntimeError:
        # already set; ignore
        pass

    if not parallel:
        job_nums = 1  # Force serial execution for debugging

    data_path = os.path.join(
        pnt.get_output_dir(), "iono_delay", "labels_csv", f"norbit_{n_orbit}_dt_{dt}s_dtrt_{dtrt}s"
    )
    save_path = os.path.join(
        pnt.get_output_dir(),
        "iono_delay",
        "raytrace_csv",
        f"norbit_{n_orbit}_dt_{dt}s_dtrt_{dtrt}s",
    )

    # Prepare argument tuples for each job
    for sat_num in sat_nums:
        for gnss_const in gnss_consts:
            for freq_family in freq_families:
                for Rz12 in Rz12s:
                    for straight_ray in straight_rays:

                        if parallel:
                            job_args = [
                                (
                                    job_id,
                                    job_nums,
                                    epoch_str,
                                    sat_num,
                                    gnss_const,
                                    freq_family,
                                    Rz12,
                                    kp,
                                    correction,
                                    straight_ray,
                                    data_path,
                                    save_path,
                                )
                                for job_id in range(job_nums)
                            ]
                            all_ok = True
                            with mp.Pool(processes=job_nums) as pool:
                                results = []
                                for ja in job_args:
                                    results.append(pool.apply_async(run_single_job, ja))

                                pool.close()

                                # Wait for all to complete and check for errors
                                for res in results:
                                    try:
                                        res.get()
                                    except Exception as e:
                                        all_ok = False
                                        print(f"[launcher] One job failed: {e}", flush=True)

                                pool.join()

                                if all_ok:
                                    print("All jobs completed successfully")
                                else:
                                    print("Some jobs FAILED — check logs above.")
                        else:
                            # Run serially for debugging
                            for job_id in range(job_nums):
                                run_single_job(
                                    job_id,
                                    job_nums,
                                    epoch_str,
                                    sat_num,
                                    gnss_const,
                                    freq_family,
                                    Rz12,
                                    kp,
                                    correction,
                                    straight_ray,
                                    data_path,
                                    save_path,
                                )
    print("All done.")
