import os


def count_lines_in_file(file_path):
    try:
        with open(file_path, "r") as f:
            return sum(1 for _ in f)
    except:
        return 0


def count_files_and_lines(root_directory):
    total_files = 0
    total_lines = 0

    for dirpath, dirnames, filenames in os.walk(root_directory):
        for filename in filenames:
            total_files += 1
            file_path = os.path.join(dirpath, filename)
            total_lines += count_lines_in_file(file_path)

    return total_files, total_lines


if __name__ == "__main__":
    directories_to_analyze = [
        os.path.join("lupnt", x)
        for x in ["core", "agents", "gnss", "dynamics", "physics", "math"]
    ]

    for root_directory in directories_to_analyze:
        total_files, total_lines = count_files_and_lines(root_directory)
        print(f"In directory '{root_directory}':")
        print(f"  Total number of files: {total_files}")
        print(f"  Total number of lines: {total_lines}")
