import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Set, Dict
import pylupnt as pnt
from itertools import count
from enum import Enum

# MARKOV DECISION PROCESS FORMULATION

# A Markov Decision Process is a general framework for modeling and solving decision making problems. The agent (the satellite) chooses an action in the current state, receives a reward, and then transitions to the next state. The process is repeated over the planning horizon. The state space S, action space A, transition function T , and reward function R, deﬁne an MDP
#       M = (S, A, T, R)
# Relevant to the Earth observing satellite tasking problem is a set of locations I, which have been requested to be collected over a time interval [0, H] with H as the planning horizon. Each image i ∈ I is deﬁned by an Earth-ﬁxed center point and has an associated reward for collection ri. Over the horizon, each image has a set of windows of opportunity Oi. Each individual opportunity for image i is denoted oi ∈ Oi and has a start time ts and end time te during which the image could be collected. The problem is to decide which images should be collected to maximize the total reward over the planning horizon.


@dataclass(frozen=True, order=True)
class Location:
    # Location on the body surface
    lat: float  # [deg] Latitude
    lon: float  # [deg] Longitude
    alt: float  # [km] Altitude
    r: np.ndarray[float, 3]  # [km] Position in body frame

    def __post_init__(self):
        self.r = pnt.geographical_to_cartesian(
            [self.lat, self.lon, self.alt],
            pnt.R_MOON,
        )


@dataclass(order=True)
class Task:
    time_start: float  # [s] Start time
    time_end: float  # [s] End time
    duration: float  # [s] Duration
    power: float  # [W] Power required
    data: float  # [Mbps] Data generated
    reward: float  # [-] Reward for task

    id: int = field(compare=False, init=False)
    id_counter: int = 0

    def __post_init__(self):
        self.id = self.id_counter
        self.id_counter += 1

    def __hash__(self):
        return hash(self.id)


# LRO max slew rate in Observing mode is 0.1 deg/sec and HGA max rate is 0.5 deg/sec).
@dataclass(frozen=True, order=True)
class Satellite:
    power_capacity: float = 100  # [Wh]
    data_capacity: float = 100  # [Mb]
    slew_power: float = 1  # [W]
    slew_rate: float = 0.25  # [deg/s]


@dataclass
class SatelliteState:
    time: float  # [s] Start time of the last action taken
    opportunities: set[Opportunity]  # [-] Future opportunities
    data: float  # [GB] Current on-board data usage
    power: float  # [Wh] Current spacecraft power


@dataclass
class SatelliteAction:
    opportunity: Opportunity  # [-] Next opportunity to take


# For each state, the action space A(s) is the set of all possible actions with start-time after the current time that are feasible for the spacecraft to transition to given the last collection or contact time tsp:
#       A(s) = {a | ts > t, C(tsp, ts) = 1}
# C(tsp, ts) is the agility constraint function that indicates whether a slew between the pointing conﬁgurations for ac- tions with start times tsp and ts is feasible. It is a binary function C(tsp, ts) ∈ {0, 1} where 1 denotes feasibility and 0 infeasibility of the transition.


def slew_time(a: SatelliteAction) -> float:
    return 0


@dataclass
class SatelliteTaskingMdp:

    # satellite: Satellite

    def transition(self, s: SatelliteState, a: SatelliteAction) -> SatelliteState:
        time = min(s.time, a.opportunity.time_start) + a.opportunity.duration
        if a.opportunity.type not in [TaskType.SUN_POINTING, TaskType.DOWNLINK]:
            opportunities = s.opportunities - {a.opportunity}
        data = s.data + a.opportunity.data * a.opportunity.duration
        power = s.power - a.opportunity.power * a.opportunity.duration
        sp = SatelliteState(time, opportunities, data, power)
        return sp

    def reward(self, s: SatelliteState, a: SatelliteAction) -> float:
        return a.opportunity.reward

    def available_actions(self, s: SatelliteState) -> List[SatelliteAction]:
        return [
            SatelliteAction(opp)
            for opp in s.opportunities
            if s.time <= opp.time_end - opp.duration
            # and power
            # and data
        ]

    def is_terminal(self, s: SatelliteState) -> bool:
        # Only sun pointing and downlink opportunities left
        return len(self.available_actions(s)) == 2


@dataclass
class SmdpForwardSearch:
    model: SatelliteTaskingMdp

    def select_actions(self, s, d, gamma):
        if d == 0 or self.model.is_terminal(s):
            return (None, 0)
        a_best, v_best = None, float("-inf")
        for action in self.model.available_actions(s):
            sp = self.model.transition(s, action)
            value = (
                self.model.reward(s, action)
                + gamma * self.select_actions(sp, d - 1, gamma)[1]
            )
            if value > v_best:
                a_best, v_best = action, value
        return (a_best, v_best)

    def solve(self, s, solve_depth, discount_factor):
        policy = []
        a, v = self.select_actions(s, solve_depth, discount_factor)
        policy.append((s, a))
        s = self.model.transition(s, a)
        while not self.model.is_terminal(s):
            a, v = self.select_actions(s, solve_depth, discount_factor)
            policy.append((s, a))
            s = self.model.transition(s, a)
        return policy
