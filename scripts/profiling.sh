#!/bin/bash

# Check if the binary is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <path_to_binary>"
  exit 1
fi

BINARY=$1
OUTPUT_DIR="profile_results"
PROF_FILE="$OUTPUT_DIR/profiler.prof"
PDF_OUTPUT="$OUTPUT_DIR/profile_report.pdf"
TEXT_OUTPUT="$OUTPUT_DIR/profile_report.txt"

# Create output directory
mkdir -p $OUTPUT_DIR

# Find libprofiler.dylib in common locations
if [ -f "/usr/local/lib/libprofiler.dylib" ]; then
  LIBPROFILER_PATH="/usr/local/lib/libprofiler.dylib"
elif [ -f "/opt/homebrew/lib/libprofiler.dylib" ]; then
  LIBPROFILER_PATH="/opt/homebrew/lib/libprofiler.dylib"
else
  echo "libprofiler.dylib not found. Please install gperftools."
  exit 1
fi

# Set up the library path
export DYLD_FALLBACK_LIBRARY_PATH="/usr/local/lib:/opt/homebrew/lib:/usr/lib:/lib:/opt/local/lib"

# Run the binary with profiling enabled
echo "Running the binary with gperftools profiling..."
export DYLD_INSERT_LIBRARIES="$LIBPROFILER_PATH"
export CPUPROFILE=$PROF_FILE
$BINARY

# Unset the environment variables
unset DYLD_INSERT_LIBRARIES
unset CPUPROFILE

# Check if the profile file was generated
if [ ! -f "$PROF_FILE" ]; then
  echo "Failed to generate profile file."
  exit 1
fi

# Generate text report using pprof
echo "Generating text report using pprof..."
pprof --text $BINARY $PROF_FILE > $TEXT_OUTPUT

# Verify the text report was generated
if [ ! -f "$TEXT_OUTPUT" ]; then
  echo "Failed to generate text report."
  exit 1
fi

# Create a PDF from the report using pandoc
echo "Creating PDF report..."
pandoc $TEXT_OUTPUT -o $PDF_OUTPUT

# Open the PDF report
if [ -f "$PDF_OUTPUT" ]; then
  echo "Opening the PDF report..."
  open "$PDF_OUTPUT"
else
  echo "Failed to generate PDF report."
fi

echo "Profiling complete. Report saved to $PDF_OUTPUT"
