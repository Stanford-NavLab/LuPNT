import os


def count_files_lines_in_dir(directory, extensions):
    """
    Recursively counts the files and lines in files with specified extensions within the given directory.

    Args:
    directory (str): The path to the directory to search.
    extensions (set): A set of file extensions to include in the count.

    Returns:
    int, int: The total number of lines and files in the directory.
    """
    total_lines = 0
    total_files = 0
    subdirectory_counts = {}

    # Walk through all subdirectories and files
    for root, dirs, files in os.walk(directory):
        # Count lines and files in the current directory
        dir_lines = 0
        dir_files = 0
        for file in files:
            if any(file.endswith(ext) for ext in extensions):
                file_path = os.path.join(root, file)
                with open(file_path, "r", encoding="utf-8") as f:
                    line_count = sum(1 for line in f)
                    dir_lines += line_count
                dir_files += 1

        # Accumulate line and file counts from subdirectories
        if root != directory:
            parent_dir = os.path.dirname(root)
            subdirectory_counts.setdefault(parent_dir, (0, 0))
            subdirectory_counts[parent_dir] = (
                subdirectory_counts[parent_dir][0] + dir_lines,
                subdirectory_counts[parent_dir][1] + dir_files,
            )

        if dir_lines > 0 or dir_files > 0:
            total_lines += dir_lines
            total_files += dir_files
            subdirectory_counts[root] = (dir_lines, dir_files)

    return total_lines, total_files, subdirectory_counts


def print_directory_counts(directory, counts, prefix=""):
    """
    Prints the file and line counts in a hierarchical structure.

    Args:
    directory (str): The path to the directory.
    counts (dict): A dictionary with paths as keys and tuples (line counts, file counts) as values.
    prefix (str): A prefix for indentation to reflect hierarchy.
    """
    lines, files = counts[directory]
    print(f"{prefix}{os.path.basename(directory)}: {files} files, {lines} lines")
    for root, dirs, files in os.walk(directory):
        for dir in sorted(dirs):
            path = os.path.join(root, dir)
            if path in counts and path != directory:
                print_directory_counts(path, counts, prefix + "    ")


# Example usage:
directory_paths = ["python", "lupnt", "examples", "scripts", "test"]
extensions = {".py", ".pyi", ".txt", ".cc", ".h", ".sh", ".cmake"}
for directory_path in directory_paths:
    total_lines, total_files, counts = count_files_lines_in_dir(
        directory_path, extensions
    )
    print_directory_counts(directory_path, counts)
    print("\n")
