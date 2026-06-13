#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
cd "${repo_root}"

config="projects/GNSS_Filtering/gnss_filtering_config.yaml"
build_dir="build-gnss-filtering"
run_delays=1
plot=1
overwrite_delays=0
serial_delays=0
delay_workers=""

usage() {
  cat <<'EOF'
Usage: projects/GNSS_Filtering/run_pipeline.sh [options]

Build and run the GNSS filtering pipeline:
  1. C++ link precompute
  2. Python GCPM delay precompute
  3. C++ Monte Carlo filtering
  4. Python post-processing plots

Options:
  --config PATH       YAML config path
  --skip-delays       Do not run the GCPM delay precompute stage. Use only when
                      the delay table already exists or plasma truth is disabled.
  --workers N         Number of Python delay worker processes
  --serial-delays     Run delay precompute in serial
  --overwrite-delays  Recompute existing delay rows
  --no-plot           Do not run post-processing plots
  -h, --help          Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      config="$2"
      shift 2
      ;;
    --skip-delays)
      run_delays=0
      shift
      ;;
    --workers)
      delay_workers="$2"
      shift 2
      ;;
    --serial-delays)
      serial_delays=1
      shift
      ;;
    --overwrite-delays)
      overwrite_delays=1
      shift
      ;;
    --no-plot)
      plot=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

sim_bin="${build_dir}/examples/ex_lunar_gnss_odts"
delay_args=("${config}")
if [[ -n "${delay_workers}" ]]; then
  delay_args+=("--workers" "${delay_workers}")
fi
if [[ "${serial_delays}" -eq 1 ]]; then
  delay_args+=("--serial")
fi
if [[ "${overwrite_delays}" -eq 1 ]]; then
  delay_args+=("--overwrite")
fi

echo "[1/5] Building GNSS filtering executables"
cmake -B "${build_dir}" -GNinja -DCMAKE_BUILD_TYPE=Release
cmake --build "${build_dir}" --target ex_lunar_gnss_odts

echo "[2/5] Precomputing GNSS links"
"${sim_bin}" --config "${config}" --precompute

if [[ "${run_delays}" -eq 1 ]]; then
  echo "[3/5] Precomputing GCPM plasma/ionosphere delays"
  python projects/GNSS_Filtering/precompute_delays.py "${delay_args[@]}"
else
  echo "[3/5] Skipping GCPM plasma/ionosphere delay precompute"
fi

echo "[4/5] Running Monte Carlo GNSS filtering"
"${sim_bin}" --config "${config}" --run

if [[ "${plot}" -eq 1 ]]; then
  echo "[5/5] Generating post-processing plots"
  python projects/GNSS_Filtering/plot_gnss_filtering_results.py
else
  echo "[5/5] Skipping post-processing plots"
fi

echo "GNSS filtering pipeline complete"
