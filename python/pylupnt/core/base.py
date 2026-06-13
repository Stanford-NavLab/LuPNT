import time
from pathlib import Path
import yaml

# BASEDIR/python/pylupnt/utils/utils.py
BASEDIR = Path(__file__).parent.parent.parent.parent


# YAML Dumper to represent lists in flow style (e.g., [a, b, c])
class BracketedListDumper(yaml.SafeDumper):
    def represent_list(self, data):
        return self.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)


BracketedListDumper.add_representer(list, BracketedListDumper.represent_list)


def enum_to_list(enum):
    enum_list = []
    for name in dir(enum):
        if not name.startswith("_"):
            enum_list.append(getattr(enum, name))
    return enum_list


def enum_to_dict(enum):
    enum_dict = {}
    for name in dir(enum):
        if not name.startswith("_"):
            enum_dict[name] = getattr(enum, name)
    return enum_dict


def get_hash(*args) -> str:
    """Generate a hash string from the given arguments."""
    import hashlib

    hash_input = "_".join(str(arg) for arg in args).encode("utf-8")
    return hashlib.md5(hash_input).hexdigest()[:8]


def get_timestamp() -> str:
    """Get current timestamp in YYYY-MM-DD HH:MM:SS format."""
    return time.strftime("%Y-%m-%d %H:%M:%S")


def convert_title_case(s: str) -> str:
    """Convert snake_case string to Title Case."""
    return " ".join(word.capitalize() for word in s.split("_"))


def set_seed(seed: int) -> None:
    import numpy as np

    """Set random seed for torch, numpy, and open3d for reproducibility."""
    np.random.seed(seed)


def wait_for_key(message: str = "Press any key to continue...", name: str = "Utils") -> None:
    """Wait for a key press (not just Enter)."""
    import pylupnt

    pylupnt.Logger.info(message, name=name)
    try:
        import msvcrt

        msvcrt.getch()
    except ImportError:
        import sys
        import tty
        import termios

        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    print()
