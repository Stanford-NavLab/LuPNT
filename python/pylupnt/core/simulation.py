import heapq
import os
import threading
import pylupnt as pnt

from datetime import datetime, timedelta
from pathlib import Path
from typing import Callable, Dict, List


from .base import BASEDIR
from pylupnt import Logger
from pylupnt import gregorian_to_time


class Event:
    def __init__(self, time: float, callback: Callable = None, priority: int = 0):
        self.time = time  # [s]
        self.callback = callback
        self.priority = priority
        self.cancelled = False

    def __lt__(self, other):
        # Sort by time then priority
        if self.time == other.time:
            return self.priority > other.priority
        return self.time < other.time


class Simulation:
    _local = threading.local()

    def __init__(self, config: pnt.Config):
        from pylupnt.agents.agent import Agent

        self.config = config
        self.name = self.__class__.__name__
        Logger.log("Setting up simulation", self.name)

        # Time
        self.time = 0.0  # [s]
        self.time_start = 0.0  # [s]
        self.duration = self.config.duration  # [s]
        if self.duration is None:
            raise ValueError("duration must be set")

        # Epoch
        if self.config.epoch_start_tai_str is not None:
            self.epoch_start_tai_str = self.config.epoch_start_tai_str
            self.epoch_start_tai = gregorian_to_time(self.epoch_start_tai_str)
            self.epoch_end_tai = (
                self.epoch_start_tai + timedelta(seconds=self.duration).total_seconds()
            )
        else:
            self.epoch_start_tai_str = None
            self.epoch_start_tai = None
            self.epoch_end_tai = None

        # Environments
        self.environments = {}

        # Agents (set it before channels!)
        self.agents: Dict[str, Agent] = {}
        if self.config.agents is not None:
            for k, v in self.config.agents.items():
                v["name"] = k
                self.agents[k] = Agent.from_config(v)
                self.agents[k].simulation = self
        else:
            self.agents = {}

        # Channels
        self.channels = {}

        # Events
        self.queue: List[Event] = []

        # Cesium
        from pylupnt.interfaces.cesium import CesiumViewer

        self.cesium_viewer = CesiumViewer()

        # Output directory
        if not os.path.isabs(self.config.output.dir):
            self.output_dir = BASEDIR / self.config.output.dir
        else:
            self.output_dir = Path(self.config.output.dir)
        if self.config.output.add_timestamp:
            self.output_dir /= datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Setup agents
        for agent in self.agents.values():
            agent.setup()

        Simulation._local.currentulation = self

    def reset(self):
        self.time = 0
        for agent in self.agents.values():
            agent.reset()
        for environment in self.environment.values():
            environment.reset()
        for channel in self.channels.values():
            channel.reset()

    @classmethod
    def from_config(cls, config: pnt.Config):
        config = pnt.load_config(config)
        return cls(config)

    def schedule(self, time: float, callback: Callable, priority: int = 0):
        if time < self.time:
            raise ValueError(f"Time must be greater than current time, got {time}")
        if time <= self.duration:
            heapq.heappush(self.queue, Event(time=time, callback=callback, priority=priority))

    def run(self):
        self.time = 0.0
        Logger.log("Simulation started", self.name, self.time)

        with Logger.tqdm(total=100, name=self.name, percentage_only=True) as pbar:
            while self.queue:
                event = heapq.heappop(self.queue)
                if event.time > self.duration:
                    Logger.warning(
                        f"Event {event.time} is greater than time_end {self.duration}",
                        self.name,
                    )
                    continue
                self.time = event.time
                event.callback(self.time)
                perc = round(self.time / self.duration * 100, 1)
                pbar.update(perc - pbar.n)

            if self.time >= self.duration:
                pbar.update(100 - pbar.n)

        Logger.log(f"Simulation finished", self.name, self.time)
