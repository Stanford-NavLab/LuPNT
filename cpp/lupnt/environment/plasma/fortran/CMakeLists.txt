# Fortran static libraries for the plasma (pecsim) module. Sources are vendored in-tree — no network
# download required.

set(FORTRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# ---- IRI2007 ----
# Override igrf.for and irifun.for with the modified versions
file(COPY ${FORTRAN_DIR}/fortran_modified/00_iri2007/igrf.for
          ${FORTRAN_DIR}/fortran_modified/00_iri2007/irifun.for
     DESTINATION ${FORTRAN_DIR}/00_iri2007
)

file(GLOB IRI2007_SOURCES "${FORTRAN_DIR}/00_iri2007/*.for")
list(FILTER IRI2007_SOURCES EXCLUDE REGEX ".*iritest\\.for$")

if(IRI2007_SOURCES)
  # Must be SHARED with hidden visibility to avoid duplicate Fortran symbols from igrf.for /
  # iridreg.for that also appear in iri2020_lib (STATIC).
  add_library(iri2007_lib SHARED ${IRI2007_SOURCES})
  target_compile_options(iri2007_lib PRIVATE -fvisibility=hidden)
else()
  message(FATAL_ERROR "IRI2007 library sources not found.")
endif()

# ---- IRI2020 ----
# Override files with modified versions
file(GLOB IRI2020_MODIFIED "${FORTRAN_DIR}/fortran_modified/00_iri2020/*.for")
foreach(f ${IRI2020_MODIFIED})
  file(COPY ${f} DESTINATION ${FORTRAN_DIR}/00_iri2020)
endforeach()

file(GLOB IRI2020_SOURCES "${FORTRAN_DIR}/00_iri2020/*.for")
list(FILTER IRI2020_SOURCES EXCLUDE REGEX ".*(iritest|irirtam-test)\\.for$")

if(IRI2020_SOURCES)
  add_library(iri2020_lib STATIC ${IRI2020_SOURCES})
else()
  message(FATAL_ERROR "IRI2020 library sources not found.")
endif()

# Copy IRI2020 data files (.asc, .dat) to the plasma data directory. Uses LUPNT_DATA_DIR (set by
# FetchLuPNTData.cmake) rather than CMAKE_SOURCE_DIR, since CMAKE_SOURCE_DIR resolves to cpp/ (not
# the repo root) when this is configured via `cmake -S cpp -B build` (e.g. `pixi run build`).
set(IRI_DATA_DEST "${LUPNT_DATA_DIR}/plasma/data/iri")
file(MAKE_DIRECTORY ${IRI_DATA_DEST})
file(GLOB IRI2020_DATA "${FORTRAN_DIR}/00_iri2020/*.asc" "${FORTRAN_DIR}/00_iri2020/*.dat")
foreach(f ${IRI2020_DATA})
  file(COPY ${f} DESTINATION ${IRI_DATA_DEST})
endforeach()

# ---- GCPM v2.4 ----
file(GLOB_RECURSE GCPM_RAW "${FORTRAN_DIR}/gcpm_v24/*.for")
list(FILTER GCPM_RAW EXCLUDE REGEX ".*(test|test2D).*\\.for$")

# Preprocess with sed: replace 'type' → 'print', remove lines starting with 'd'
set(PREPROCESSED_DIR "${CMAKE_CURRENT_BINARY_DIR}/preprocessed")
file(MAKE_DIRECTORY ${PREPROCESSED_DIR})

set(GCPM_PROCESSED)
foreach(src ${GCPM_RAW})
  get_filename_component(fname ${src} NAME)
  set(out "${PREPROCESSED_DIR}/${fname}")
  add_custom_command(
    OUTPUT ${out}
    COMMAND sed -e "s/type/print/g" -e "/^d/ d" ${src} > ${out}
    DEPENDS ${src}
    COMMENT "Preprocessing ${fname}"
  )
  list(APPEND GCPM_PROCESSED ${out})
endforeach()

if(NOT GCPM_PROCESSED)
  message(FATAL_ERROR "GCPM library sources not found.")
endif()

set_source_files_properties(
  ${GCPM_PROCESSED} PROPERTIES COMPILE_FLAGS "-ffixed-form -std=legacy -fno-automatic"
)
add_library(gcpm_lib STATIC ${GCPM_PROCESSED})

# ---- xform ----
file(GLOB_RECURSE XFORM_SOURCES "${FORTRAN_DIR}/xform/*.for" "${FORTRAN_DIR}/xform/*.FOR")

if(XFORM_SOURCES)
  add_library(xform_lib STATIC ${XFORM_SOURCES})
else()
  message(FATAL_ERROR "xform library sources not found.")
endif()

# ---- IGRF14 ----
file(GLOB IGRF_SOURCES "${FORTRAN_DIR}/igrf14/igrf14.f")

if(IGRF_SOURCES)
  add_library(igrf_lib STATIC ${IGRF_SOURCES})
else()
  message(FATAL_ERROR "IGRF14 library sources not found.")
endif()
