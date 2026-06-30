# Downloads and extracts the LuPNT_data reference-data archive (ephemeris, GNSS antenna/clock
# products, plasma coefficients, TLEs) into <repo>/data/LuPNT_data if it isn't already present.
# The archive is gitignored (see .gitignore) and not checked into the repo; CI fetches it from the
# same source (see .github/workflows/{ubuntu,macos,python,examples,install}.yml).
#
# Anchored via CMAKE_CURRENT_LIST_DIR (this file's own location) rather than CMAKE_SOURCE_DIR,
# since CMAKE_SOURCE_DIR differs depending on which CMakeLists.txt is the top-level project for a
# given `pixi run` task (e.g. `cmake -S cpp -B build` has CMAKE_SOURCE_DIR=cpp/, not the repo
# root).

get_filename_component(LUPNT_REPO_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(LUPNT_DATA_DIR "${LUPNT_REPO_ROOT}/data/LuPNT_data")

option(LUPNT_FETCH_DATA
       "Automatically download data/LuPNT_data (ephemeris/GNSS/plasma/TLE reference data) if missing"
       ON
)

if(LUPNT_FETCH_DATA)
  set(LUPNT_DATA_MARKER "${LUPNT_DATA_DIR}/gnss/igs20.atx")

  if(NOT EXISTS "${LUPNT_DATA_MARKER}")
    message(
      STATUS
        "LuPNT data: ${LUPNT_DATA_MARKER} not found; downloading data/LuPNT_data (~650MB, one-time)..."
    )

    set(LUPNT_DATA_ZIP "${CMAKE_BINARY_DIR}/LuPNT_data.zip")
    file(DOWNLOAD "https://bit.ly/LuPNT_data" "${LUPNT_DATA_ZIP}" SHOW_PROGRESS
         STATUS LUPNT_DATA_DOWNLOAD_STATUS
    )
    list(GET LUPNT_DATA_DOWNLOAD_STATUS 0 LUPNT_DATA_DOWNLOAD_CODE)

    if(NOT LUPNT_DATA_DOWNLOAD_CODE EQUAL 0)
      list(GET LUPNT_DATA_DOWNLOAD_STATUS 1 LUPNT_DATA_DOWNLOAD_MSG)
      file(REMOVE "${LUPNT_DATA_ZIP}")
      message(
        WARNING
          "LuPNT data: download failed (${LUPNT_DATA_DOWNLOAD_MSG}). Download it manually from "
          "https://bit.ly/LuPNT_data and extract its LuPNT_data/ folder into "
          "${LUPNT_REPO_ROOT}/data/, or configure with -DLUPNT_FETCH_DATA=OFF to silence this check."
      )
    else()
      file(REMOVE_RECURSE "${LUPNT_DATA_DIR}")
      file(MAKE_DIRECTORY "${LUPNT_REPO_ROOT}/data")
      file(ARCHIVE_EXTRACT INPUT "${LUPNT_DATA_ZIP}" DESTINATION "${LUPNT_REPO_ROOT}/data")
      file(REMOVE "${LUPNT_DATA_ZIP}")
      message(STATUS "LuPNT data: extracted to ${LUPNT_DATA_DIR}")
    endif()
  endif()
endif()
