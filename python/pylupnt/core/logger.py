import re
import shutil
import sys
import time

from tqdm import tqdm
from typing import Union


# ANSI color codes
class Colors:
    RESET = "\033[0m"
    CYAN = "\033[36m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    RED = "\033[31m"
    BLUE = "\033[34m"
    WHITE = "\033[37m"


def format_time(seconds: float, short: bool = False) -> str:
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = seconds % 60
    if not short:
        return f"{hours:02d}:{minutes:02d}:{secs:05.2f}"

    msg = ""
    msg += f"{hours:02d}:" if hours > 0 else ""
    msg += f"{minutes:02d}:" if minutes > 0 else ""
    msg += f"{secs:04.1f}"
    return msg


class Commands:
    MOVE_UP = "\033[1A"
    CLEAR_LINE = "\033[K"


class Logger:
    time_start = time.time()
    time_color = Colors.CYAN
    name_color = Colors.GREEN
    pbar = None

    DEBUG = 0
    INFO = 1
    WARNING = 2
    ERROR = 3
    log_level = INFO

    log_level_map = {
        "DEBUG": 0,
        "INFO": 1,
        "WARNING": 2,
        "ERROR": 3,
    }

    @staticmethod
    def set_log_level(log_level: Union[int, str]):
        if isinstance(log_level, str):
            Logger.log_level = Logger.log_level_map[log_level]
        else:
            Logger.log_level = log_level

    @staticmethod
    def clear_pbars():
        for pbar in Logger.pbars:
            if pbar.n >= pbar.total:
                pbar.close()
        Logger.pbars = [pbar for pbar in Logger.pbars if not pbar.n >= pbar.total]

    @staticmethod
    def log(
        message: str,
        name: str = None,
        time_sim: float = None,
        style: str = Colors.WHITE,
    ):
        elapsed_real = format_time(time.time() - Logger.time_start, short=True)

        msg = f"[{Logger.time_color}{elapsed_real}{Colors.RESET}]"
        msg_len = len(elapsed_real) + 2
        if name is not None:
            msg += f"[{Logger.name_color}{name}{Colors.RESET}]"
            msg_len += len(name) + 2
        if time_sim is not None:
            elapsed_sim = format_time(time_sim, short=True)
            msg += f"[{Logger.time_color}{elapsed_sim}{Colors.RESET}]"
            msg_len += len(elapsed_sim) + 2
        if message is not None:
            msg += f" {style}{message}{Colors.RESET}"
            msg_len += len(message) + 1

        if Logger.pbar is not None:
            pad = " " * (shutil.get_terminal_size().columns - msg_len)
            # Update the pbar desc with the new elapsed real
            new_desc = f"[{Logger.time_color}{elapsed_real}{Colors.RESET}]"
            Logger.pbar.desc = re.sub(r"^\[.*?\]", new_desc, Logger.pbar.desc)
            Logger.pbar.write(msg + pad)
            Logger.pbar.refresh()
        else:
            print(msg)

    @staticmethod
    def tqdm(sequence=None, **kwargs):
        elapsed_real = format_time(time.time() - Logger.time_start, short=True)
        percentage_only = kwargs.pop("percentage_only", False)
        if not percentage_only:
            l_bar = "{desc} {percentage:3.0f}%|"
            bar = "{bar}"
            r_bar = "| {n_fmt}/{total_fmt} [{elapsed}<{remaining},{rate_fmt}{postfix}]"
            bar_format = f"{l_bar}{bar}{r_bar}"
        else:
            bar_format = "{desc}|{bar}|{percentage:3.0f}% [{remaining} left]"

        msg = f"[{Logger.time_color}{elapsed_real}{Colors.RESET}]"
        name = kwargs.pop("name", None)
        if name is not None:
            msg += f"[{Logger.name_color}{name}{Colors.RESET}]"
        desc = kwargs.pop("desc", None)
        if desc is not None:
            msg += f" {desc}"
        if sequence is None:
            pbar = tqdm(desc=msg, bar_format=bar_format, **kwargs, file=sys.stdout)
        else:
            pbar = tqdm(sequence, desc=msg, bar_format=bar_format, **kwargs, file=sys.stdout)

        Logger.pbar = pbar
        return pbar

    @staticmethod
    def info(message: str, name: str = None, time_sim: float = None):
        if Logger.log_level <= Logger.INFO:
            Logger.log(message, name, time_sim, style=Colors.WHITE)

    @staticmethod
    def warning(message: str, name: str = None, time_sim: float = None):
        if Logger.log_level <= Logger.WARNING:
            Logger.log(message, name, time_sim, style=Colors.YELLOW)

    @staticmethod
    def error(message: str, name: str = None, time_sim: float = None):
        if Logger.log_level <= Logger.ERROR:
            Logger.log(message, name, time_sim, style=Colors.RED)

    @staticmethod
    def debug(message: str, name: str = None, time_sim: float = None):
        if Logger.log_level <= Logger.DEBUG:
            Logger.log(message, name, time_sim, style=Colors.BLUE)

    @staticmethod
    def reset_start_time():
        Logger.time_start = time.time()
