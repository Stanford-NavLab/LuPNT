import argparse
from pathlib import Path
from git import Repo
import sys


def parse_args():
    parser = argparse.ArgumentParser("Manage the project version including git tags.")
    parser.add_argument(
        "--bump",
        default="none",
        const="none",
        nargs="?",
        choices=["major", "minor", "patch", "none"],
        help="Bump version of project",
    )
    parser.add_argument(
        "--version_file", default="source/version.txt", help="Path to version file"
    )
    parser.add_argument(
        "--tag", dest="tag", action="store_true", default=False, help="Create git tag"
    )
    args = parser.parse_args()

    assert Path(args.version_file).exists()

    return args


def bump(bump_type: str = "patch", version_file="source/version.txt"):
    with open(version_file, "r") as f:
        lines = f.readlines()
    version = [int(line.split(" ")[1].strip()) for line in lines if line]

    # Increment the version, resetting minor/patches
    if bump_type == "major":
        version[0] += 1
        version[1] = 0
        version[2] = 0
    if bump_type == "minor":
        version[1] += 1
        version[2] = 0
    if bump_type == "patch":
        version[2] += 1

    if bump_type != "none":
        with open(version_file, "w") as f:
            f.write("MAJOR {}\nMINOR {}\nPATCH {}\n".format(*version))

        # Update pyproject.toml
        update_pyproject_toml(version)
    else:
        pass
        # print("Not incrementing version, bump was 'none'.\n")

    version_str = ".".join(map(str, version))
    return version_str


def update_pyproject_toml(version):
    pyproject_path = Path(__file__).resolve().parent.parent / "pyproject.toml"
    if not pyproject_path.exists():
        print("pyproject.toml not found.")
        return

    with open(pyproject_path, "r") as f:
        lines = f.readlines()

    with open(pyproject_path, "w") as f:
        for line in lines:
            if line.startswith("version ="):
                f.write(f'version = "{version[0]}.{version[1]}.{version[2]}"\n')
            else:
                f.write(line)


def create_tag(version_str, message="Tagged for release"):
    version_str = "v{}".format(version_str)
    repo = Repo(".")
    if repo.is_dirty():
        resp: str = input("Repo is dirty do you want to continue? [y/n]: ")
        if resp.lower() != "y":
            print("Quiting...")
            sys.exit()
    hexsha = repo.head.commit.hexsha
    print("Creating new tag: {} ; commit: {}".format(version_str, hexsha))
    new_tag = repo.create_tag(version_str, message=message)
    return version_str


def main():
    args = parse_args()
    if args.tag:
        version_str = bump(bump_type="none", version_file=args.version_file)
        new_tag = create_tag(version_str)
        print("Don't forget to push the tags with:")
        print("\t git push origin {}  OR".format(new_tag))
        print("\t git push origin --tags")
    else:
        version_str = bump(bump_type=args.bump, version_file=args.version_file)
        print("Incremented version to {}".format(version_str))


if __name__ == "__main__":
    main()
