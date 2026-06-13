#!/usr/bin/env python3
"""Check formatting for files changed from a Git base ref."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


CPP_EXTENSIONS = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx"}
CMAKE_EXTENSIONS = {".cmake"}


def run(args: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, check=check, text=True, capture_output=True)


def git_lines(args: list[str]) -> list[str]:
    result = run(["git", *args])
    return [line for line in result.stdout.splitlines() if line]


def changed_files(base_ref: str) -> list[str]:
    diff_files = git_lines(["diff", "--name-only", "--diff-filter=ACMRT", base_ref, "--"])
    untracked_files = git_lines(["ls-files", "--others", "--exclude-standard"])
    return sorted(set(diff_files + untracked_files))


def is_cmake_file(path: str) -> bool:
    name = Path(path).name
    return name == "CMakeLists.txt" or Path(path).suffix in CMAKE_EXTENSIONS


def main() -> int:
    base_ref = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("FORMAT_BASE_REF", "HEAD")
    files = changed_files(base_ref)
    cpp_files = [path for path in files if Path(path).suffix in CPP_EXTENSIONS]
    cmake_files = [path for path in files if is_cmake_file(path)]

    failed = False
    if cpp_files:
        result = run(["clang-format", "--dry-run", "--Werror", *cpp_files], check=False)
        if result.returncode:
            sys.stderr.write(result.stderr)
            failed = True

    if cmake_files:
        result = run(["cmake-format", "--check", *cmake_files], check=False)
        if result.returncode:
            sys.stdout.write(result.stdout)
            sys.stderr.write(result.stderr)
            failed = True

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
