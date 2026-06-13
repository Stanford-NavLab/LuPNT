import numpy as np
import pylupnt as pnt
from src.constellation_design import create_x, N_PHASE_VEC, int_to_sma


def create_initial_population(problem, sat_range_phases, pop_size, use_float_sma, seed=1):
    """
    Create an initial population for the optimization problem.

    Parameters:
    problem : pymoo.core.problem.Problem
        The optimization problem instance.
    pop_size : int
        The desired population size.

    Returns:
    X0 : np.ndarray
        An array of shape (pop_size, n_var) representing the initial population.
    """
    np.random.seed(seed)

    x = []
    num_x = 0

    family_ratios = [1 / 4, 1 / 4, 1 / 2]  # ratio of each family: North+South, Circular, Hybrid
    printed_len = []

    while len(x) < pop_size:
        if np.remainder(num_x, 10) == 0 and num_x not in printed_len:
            print("Creating individual {}/{}".format(num_x + 1, pop_size))
            printed_len.append(num_x)

        # set new deployment satellite numbers for each satellite phase
        phase1_sat = np.random.randint(sat_range_phases[0][0], sat_range_phases[0][1] + 1)
        phase2_sat = np.random.randint(
            max(sat_range_phases[1][0] - phase1_sat, 0), sat_range_phases[1][1] + 1 - phase1_sat
        )
        phase3_sat = np.random.randint(
            max(sat_range_phases[2][0] - phase1_sat - phase2_sat, 0),
            sat_range_phases[2][1] + 1 - phase1_sat - phase2_sat,
        )
        tot_sat_used = phase1_sat + phase2_sat + phase3_sat

        # decide which family to use
        if num_x < (family_ratios[0] * pop_size):
            family = 0
        elif num_x < ((sum(family_ratios[:2]) * pop_size)):
            family = 1
        else:
            family = 2

        # Design ELFO and CLFOs --------------------------------------------------------------------------------------
        # sma
        if use_float_sma:
            a_elfo = np.random.uniform(6000, 15000)  # a elfo ELFO for pole coverage
            a_clfo = np.random.uniform(6000, 15000)  # a clfo ELFO for equatorial coverage

            a_elfo_km = a_elfo
        else:
            a_elfo = np.random.randint(0, 6)  # a random integer between 0 and 6
            a_clfo = np.random.randint(0, 6)  # a random integer between 0 and 6

            a_elfo_km = int_to_sma(a_elfo)

        # ecc
        ecc_elfo_max = min(0.7, 1 - (pnt.R_MOON + 100) / a_elfo_km)
        ecc_clfo = 0.001  # almost circular
        ecc_elfo = np.random.uniform(0.3, ecc_elfo_max)

        # Plane and number of satellites per plane ------------------------------
        # ELFO South
        if family == 0 or family == 2:  # North+South or Hybrid
            nplane_elfo_s = np.random.randint(2, 6)
            min_sat_plane = int(np.ceil(tot_sat_used / nplane_elfo_s))
            nsat_plane_elfo_s = np.random.randint(
                max(min_sat_plane, 2), max(8, min_sat_plane + 1)
            )  # make sure walker 1 is elfo to fit all phase 1 sats
        else:
            nplane_elfo_s = 2  # minimum 2 planes (not used)
            nsat_plane_elfo_s = 2  # minimum 2 sats per plane

        # CLFO
        if family == 1:
            nplane_clfo = np.random.randint(3, 8)
            min_sat_plane = int(np.ceil(tot_sat_used / nplane_clfo))
            # print("clfo  plane:{}  sat/plane:{}   total sats:{}".format(nplane_clfo, min_sat_plane, tot_sat_used))
            nsat_plane_clfo = np.random.randint(
                max(min_sat_plane, 2), max(8, min_sat_plane + 1)
            )  # make sure walker 2 is clfo to fit all satellies
        elif family == 0:
            nplane_clfo = 2  # minimum 2 planes (not used)
            nsat_plane_clfo = 2  # minimum 2 sats per plane
        else:
            nplane_clfo = np.random.randint(2, 6)
            nsat_plane_clfo = np.random.randint(2, 8)

        # ELFO North
        if family == 2:  # North+South or Hybrid
            nplane_elfo_n = np.random.randint(2, 6)
            sat_walker1 = int(nplane_elfo_s * nsat_plane_elfo_s)
            sat_walker3_min = tot_sat_used - sat_walker1
            min_sat_plane = int(np.ceil(sat_walker3_min / nplane_elfo_n))
            nsat_plane_elfo_n = np.random.randint(
                max(min_sat_plane, 2), max(8, min_sat_plane + 1)
            )  # make sure walker 3 is elfo to fit all sats
        elif family == 1:  # Circular
            nplane_elfo_n = 2  # minimum 2 planes (not used)
            nsat_plane_elfo_n = 2
        else:  # North+South (= use symmetric configuration)
            nplane_elfo_n = nplane_elfo_s
            nsat_plane_elfo_n = nsat_plane_elfo_s

        # phasing angles -------------------------------------------------
        phasing_elfo_s = np.random.randint(1, 2)
        phasing_clfo = np.random.randint(1, 2)
        phasing_elfo_n = np.random.randint(1, 2)

        Omega0 = np.random.uniform(0, 2 * np.pi)
        Omega1 = np.random.uniform(0, 2 * np.pi)
        Omega2 = np.random.uniform(0, 2 * np.pi)

        walker_params = [
            np.array(
                [a_elfo, ecc_elfo, 1, nplane_elfo_s, nsat_plane_elfo_s, phasing_elfo_s, Omega0]
            ),  # walker 0 (w=90, elfo)
            np.array(
                [a_clfo, ecc_clfo, 0, nplane_clfo, nsat_plane_clfo, phasing_clfo, Omega1]
            ),  # walker 1 (w=-90, clfo)
            np.array(
                [a_elfo, ecc_elfo, 0, nplane_elfo_n, nsat_plane_elfo_n, phasing_elfo_n, Omega2]
            ),  # walker 2 (w=-90, elfo)
        ]

        walker1_sat = int(nplane_elfo_s * nsat_plane_elfo_s)
        walker2_sat = int(nplane_clfo * nsat_plane_clfo)
        walker3_sat = int(nplane_elfo_n * nsat_plane_elfo_n)
        tot_sat_walker = walker1_sat + walker2_sat + walker3_sat

        if tot_sat_walker > N_PHASE_VEC:
            # total satellites of walker exceeds the limit
            # print("Total sats of walker {} exceeds N_PHASE_VEC".format(tot_sat_walker))
            continue

        if tot_sat_walker < tot_sat_used:
            # total satellites of walker is less than the number of sats
            # print("Total sats of walker {} is less than the number of sats {} to be used".format(tot_sat_walker, tot_sat_used))
            continue

        # phase allocation --------------------------------------------------------------
        phase_sats = np.array([phase1_sat, phase2_sat, phase3_sat])
        walker_sats = np.array([walker1_sat, walker2_sat, walker3_sat])

        if family == 0:  # North+South
            preference_walker = [0, 2]  # prefer walker 0 and walker 2
        elif family == 1:  # Circular
            preference_walker = [1]  # prefer walker 1
        elif family == 2:  # Hybrid
            preference_walker = [0, 1, 2]  # prefer walker 0 and walker 1

        x_phase, success = allocate_sats_to_stage(walker_sats, phase_sats, preference_walker)

        if not success:  # failed to allocate the sats
            continue

        x_vec = create_x(walker_params, x_phase)
        x_dict = problem._x_from_vector(x_vec)
        x.append(x_dict)
        num_x += 1

    print("Created initial population with {} individuals".format(len(x)))
    print(x[0])

    return x


def allocate_sats_to_stage(walker_sats, phase_sats, preference_walker, debug=False):
    x_phase = -1 * np.ones((N_PHASE_VEC,), dtype=object)
    total_sats = np.sum(phase_sats)
    total_slot_num = np.sum(walker_sats)
    n_phase = len(phase_sats)

    usable_slots = []
    walker_slot_idx = {}
    for w in preference_walker:
        start_idx = np.sum(walker_sats[:w])
        end_idx = start_idx + walker_sats[w]
        if debug:
            print(f"walker {w}: slots {start_idx} to {end_idx-1} (total {walker_sats[w]})")
        usable_slots.extend(list(np.arange(start_idx, end_idx)))
        walker_slot_idx[w] = list(np.arange(start_idx, end_idx))

    success = True

    if total_sats > len(usable_slots):
        # print(walker_slot_idx)
        # print(f"Not enough slots for the preferred walkers: total sats {total_sats}, usable slots {len(usable_slots)}")
        return x_phase, False

    for phase in range(3):
        phase_remain_nsat = phase_sats[phase]
        for w in preference_walker:
            # check if the walker is already full
            empty_slots_w = np.where(x_phase[walker_slot_idx[w]] == -1)[0]
            empty_slots_w_num = empty_slots_w.size
            if empty_slots_w_num > 0:  # there is still empty slot
                # allocate as many as possible to the preferred walker
                allocatable_num = min(phase_remain_nsat, empty_slots_w_num)
                if debug:
                    print(
                        f"Remaining sats for phase {phase}: {phase_remain_nsat}, empty slots in walker {w}: {empty_slots_w_num}, allocatable: {allocatable_num}"
                    )
                if allocatable_num > 0:
                    chosen_idx = np.random.choice(
                        empty_slots_w, size=allocatable_num, replace=False
                    )
                    idx = [int(walker_slot_idx[w][ci]) for ci in chosen_idx]
                    x_phase[idx] = phase
                    phase_remain_nsat -= allocatable_num

            if phase_remain_nsat == 0:  # all sats for the phase are allocated
                break

        if phase_remain_nsat > 0:
            # there is a left-over sats that could not be allocated to the preferred walkers
            success = False
            # print(f"Failed to allocate all sats for phase {phase}, remaining: {phase_remain_nsat}")
            break

    return x_phase, success
