import os
import nbformat

# Example usage:
directories = ["python", "lupnt", "examples", "scripts", "test"]
file_extensions = {".py", ".pyi", ".ipynb", ".txt", ".cc", ".h", ".sh", ".cmake"}


def count_files_and_lines(root_directory, depth=0):
    total_files = 0
    total_lines = 0
    items = os.listdir(root_directory)

    # Accumulate files and line counts for current directory
    for item in items:
        path = os.path.join(root_directory, item)
        if os.path.isdir(path):
            files, lines = count_files_and_lines(path, depth + 1)
            total_files += files
            total_lines += lines
        elif item.endswith(".ipynb") and ".ipynb" in file_extensions:
            try:
                with open(path, "r", encoding="utf-8") as file:
                    nb = nbformat.read(file, as_version=4)
                    for cell in nb.cells:
                        if cell.cell_type == "code":
                            line_count = (
                                cell.source.count("\n") + 1
                            )  # Count lines of code
                            total_lines += line_count
                    total_files += 1
            except Exception as e:
                print(f"Failed to read {path}: {str(e)}")
        elif any(item.endswith(ext) for ext in file_extensions):
            try:
                with open(path, "r") as file:
                    lines = file.readlines()
                    line_count = len(lines)
                    total_lines += line_count
                    total_files += 1
            except Exception as e:
                print(f"Failed to read {path}: {str(e)}")

    # Print the summary for the current directory if it has files
    if total_files > 0:
        indent = "    " * depth  # Control the level of indentation
        print(
            f"{indent}{os.path.basename(root_directory)}: {total_files} files, {total_lines} lines"
        )

    return total_files, total_lines


if __name__ == "__main__":
    total_files, total_lines = 0, 0
    for directory in directories:
        if os.path.exists(directory):
            files, lines = count_files_and_lines(directory)
            total_files += files
            total_lines += lines
            print("-" * 50)
        else:
            print(f"Directory does not exist: {directory}")

    print(f"Total: {total_files} files, {total_lines} lines")
