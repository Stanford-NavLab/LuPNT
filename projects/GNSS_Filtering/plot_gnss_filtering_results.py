from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _runs(output_dir: Path):
    files = sorted(output_dir.glob("trajectory_mc*.csv"))
    if not files:
        raise FileNotFoundError(f"No trajectory files found under {output_dir.resolve()}")
    runs = [pd.read_csv(path) for path in files]
    for df in runs:
        df["time_min"] = df["t"] / 60.0
    return runs


def _plot_components(runs, axes, columns, sigma_columns, ylabel):
    for ax, col, sig_col in zip(axes, columns, sigma_columns):
        for df in runs:
            label = f"MC {int(df['mc'].iloc[0])}"
            t = df["time_min"].to_numpy()
            err = df[col].to_numpy()
            sig = df[sig_col].to_numpy()
            ax.plot(t, err, linewidth=1.0, alpha=0.85, label=label)
            ax.fill_between(t, -sig, sig, alpha=0.10)
        ax.set_ylabel(ylabel)
        ax.set_title(col)
        ax.grid(True, alpha=0.35)


def main():
    output_dir = Path("output/gnss_filtering")
    runs = _runs(output_dir)

    fig, axes = plt.subplots(7, 1, figsize=(12, 15), sharex=True)
    _plot_components(
        runs,
        axes[0:3],
        ["r_error_m", "t_error_m", "n_error_m"],
        ["r_3sigma_m", "t_3sigma_m", "n_3sigma_m"],
        "position [m]",
    )
    _plot_components(
        runs,
        axes[3:6],
        ["rdot_error_mps", "tdot_error_mps", "ndot_error_mps"],
        ["rdot_3sigma_mps", "tdot_3sigma_mps", "ndot_3sigma_mps"],
        "velocity [m/s]",
    )

    for df in runs:
        axes[6].step(
            df["time_min"],
            df["num_tracked_satellites"],
            where="post",
            linewidth=1.1,
            label=f"MC {int(df['mc'].iloc[0])}",
        )
    axes[6].set_ylabel("tracked sats")
    axes[6].set_xlabel("receiver app elapsed coordinate time [min]")
    axes[6].set_title("Tracked Satellites")
    axes[6].grid(True, alpha=0.35)
    axes[0].legend(loc="best", ncols=2)

    fig.suptitle("GNSS Filtering RTN Errors and Tracking Count", y=0.995)
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / "gnss_filtering_rtn_tracking.png"
    fig.savefig(out, dpi=180)
    print(out.resolve())

    delay_path = output_dir / "precomputed_delays.csv"
    if not delay_path.exists():
        delay_path = output_dir / "precomputed_links.csv"
    if delay_path.exists():
        delays = pd.read_csv(delay_path)
        delays["time_min"] = (delays["t_tdb"] - delays["t_tdb"].min()) / 60.0
        fig_delay, ax_delay = plt.subplots(figsize=(12, 5))
        for (const, prn, freq), group in delays.groupby(["gnss_const", "prn", "frequency"]):
            label = f"{const}-{int(prn):02d} {freq}"
            ax_delay.plot(
                group["time_min"], group["ionosphere_plasma_delay_m"], linewidth=1.0, label=label
            )
        ax_delay.set_xlabel("receiver app elapsed coordinate time [min]")
        ax_delay.set_ylabel("plasma/ionosphere delay [m]")
        ax_delay.set_title("Precomputed Plasma/Ionosphere Delays")
        ax_delay.grid(True, alpha=0.35)
        ax_delay.legend(loc="best", ncols=2)
        fig_delay.tight_layout()
        delay_out = output_dir / "gnss_filtering_plasma_delays.png"
        fig_delay.savefig(delay_out, dpi=180)
        print(delay_out.resolve())


if __name__ == "__main__":
    main()
