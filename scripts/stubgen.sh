# usage: pybind11-stubgen [-h] [-o OUTPUT_DIR] [--root-suffix ROOT_SUFFIX] [--ignore-invalid-expressions REGEX]
#                         [--ignore-invalid-identifiers REGEX] [--ignore-unresolved-names REGEX] [--ignore-all-errors]
#                         [--enum-class-locations REGEX:LOC]
#                         [--numpy-array-wrap-with-annotated | --numpy-array-use-type-var | --numpy-array-remove-parameters]
#                         [--print-invalid-expressions-as-is] [--print-safe-value-reprs REGEX] [--exit-code] [--dry-run]
#                         [--stub-extension EXT]
#                         MODULE_NAME

# Generates stubs for specified modules

# positional arguments:
#   MODULE_NAME           module name

# options:
#   -h, --help            show this help message and exit
#   -o OUTPUT_DIR, --output-dir OUTPUT_DIR
#                         The root directory for output stubs
#   --root-suffix ROOT_SUFFIX
#                         Top-level module directory suffix
#   --ignore-invalid-expressions REGEX
#                         Ignore invalid expressions matching REGEX
#   --ignore-invalid-identifiers REGEX
#                         Ignore invalid identifiers matching REGEX
#   --ignore-unresolved-names REGEX
#                         Ignore unresolved names matching REGEX
#   --ignore-all-errors   Ignore all errors during module parsing
#   --enum-class-locations REGEX:LOC
#                         Locations of enum classes in <enum-class-name-regex>:<path-to-class> format. Example: `MyEnum:foo.bar.Baz`
#   --numpy-array-wrap-with-annotated
#                         Replace numpy/scipy arrays of 'ARRAY_T[TYPE, [*DIMS], *FLAGS]' format with 'Annotated[ARRAY_T, TYPE,
#                         FixedSize|DynamicSize(*DIMS), *FLAGS]'
#   --numpy-array-use-type-var
#                         Replace 'numpy.ndarray[numpy.float32[m, 1]]' with 'numpy.ndarray[tuple[M, typing.Literal[1]],
#                         numpy.dtype[numpy.float32]]'
#   --numpy-array-remove-parameters
#                         Replace 'numpy.ndarray[...]' with 'numpy.ndarray'
#   --print-invalid-expressions-as-is
#                         Suppress the replacement with '...' of invalid expressionsfound in annotations
#   --print-safe-value-reprs REGEX
#                         Override the print-safe check for values matching REGEX
#   --exit-code           On error exits with 1 and skips stub generation
#   --dry-run             Don't write stubs. Parses module and report errors
#   --stub-extension EXT  The file extension of the generated stubs. Must be 'pyi' (default) or 'py'

pybind11-stubgen pylupnt --output-dir "$(dirname "${BASH_SOURCE[0]}")/../source/python" --enum-class-locations "Frame:pylupnt"
# --print-invalid-expressions-as-is
