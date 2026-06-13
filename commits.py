import subprocess
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from tqdm import tqdm
import numpy as np

# Get all commits with dates
result = subprocess.run(
    ["git", "log", "--pretty=format:%H|%ad", "--date=short"],
    stdout=subprocess.PIPE,
    text=True,
)

commits = []
now = datetime.today()
cutoff = now - timedelta(days=30)

for line in result.stdout.splitlines():
    sha, date_str = line.split("|")
    date = datetime.strptime(date_str, "%Y-%m-%d")
    if date >= cutoff:
        commits.append((sha, date))

if not commits:
    print("Found 0 commits in last 30 days, switching to 'guillemc' branch.")
    subprocess.run(
        ["git", "checkout", "--quiet", "guillemc"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    # Re-run git log on guillemc branch
    result = subprocess.run(
        ["git", "log", "--pretty=format:%H|%ad", "--date=short"],
        stdout=subprocess.PIPE,
        text=True,
    )
    for line in result.stdout.splitlines():
        sha, date_str = line.split("|")
        date = datetime.strptime(date_str, "%Y-%m-%d")
        if date >= cutoff:
            commits.append((sha, date))
    print(f"Found {len(commits)} commits in last 30 days on 'guillemc' branch.")
else:
    print(f"Found {len(commits)} commits in last 30 days")

lines_added_data = []
lines_modified_data = []

# Save current HEAD to restore later
current_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()

try:
    rev_commits = list(reversed(commits))  # chronological order
    for i, (sha, date) in enumerate(tqdm(rev_commits, desc="Processing commits")):
        parent = rev_commits[i - 1][0] if i > 0 else None
        subprocess.run(
            ["git", "checkout", "--quiet", sha],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        try:
            if parent:
                diff_cmd = ["git", "diff", "--numstat", parent, sha]
            else:
                diff_cmd = ["git", "show", "--numstat", sha]
            diff_output = subprocess.check_output(diff_cmd, text=True, stderr=subprocess.DEVNULL)
            lines_added = 0
            lines_modified = 0
            for line in diff_output.splitlines():
                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        added = int(parts[0]) if parts[0] != "-" else 0
                        deleted = int(parts[1]) if parts[1] != "-" else 0
                        lines_added += added
                        lines_modified += min(added, deleted)
                    except ValueError:
                        continue
            lines_added_data.append((date, lines_added))
            lines_modified_data.append((date, lines_modified))
        except Exception:
            lines_added_data.append((date, 0))
            lines_modified_data.append((date, 0))
finally:
    subprocess.run(
        ["git", "checkout", "--quiet", current_sha],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def remove_outliers(data):
    if not data:
        return data
    values = np.array([v for _, v in data])
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1
    upper = q3 + 1.5 * iqr
    return [(d, v) for (d, v) in data if v <= upper]


if lines_added_data and lines_modified_data:
    lines_added_data = remove_outliers(lines_added_data)
    lines_modified_data = remove_outliers(lines_modified_data)
    if lines_added_data and lines_modified_data:
        # Only plot for dates present in both datasets
        dates_added, lines_added = zip(*lines_added_data)
        dates_modified, lines_modified = zip(*lines_modified_data)
        # Find intersection of dates
        common_dates = sorted(set(dates_added) & set(dates_modified))
        lines_added_dict = dict(lines_added_data)
        lines_modified_dict = dict(lines_modified_data)
        plot_dates = []
        plot_lines_added = []
        plot_lines_modified = []
        for d in common_dates:
            plot_dates.append(d)
            plot_lines_added.append(lines_added_dict.get(d, 0))
            plot_lines_modified.append(lines_modified_dict.get(d, 0))
        plt.figure(figsize=(10, 5))
        plt.plot(plot_dates, plot_lines_added, marker="o", linestyle="-", label="Lines Added")
        plt.plot(
            plot_dates,
            plot_lines_modified,
            marker="s",
            linestyle="-",
            label="Lines Modified",
        )
        plt.xlabel("Date")
        plt.ylabel("Lines")
        plt.title("Lines Added and Modified Over Last 30 Days (Outliers Removed)")
        plt.grid(True)
        plt.tight_layout()
        plt.legend()
        plt.show()
    else:
        print("No data to plot after removing outliers.")
else:
    print("No data to plot.")
