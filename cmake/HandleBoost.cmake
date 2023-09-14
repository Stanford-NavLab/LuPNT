# Store these in variables so they are automatically replicated in GTSAMConfig.cmake and such.
set(BOOST_FIND_MINIMUM_VERSION 1.65)
set(BOOST_FIND_MINIMUM_COMPONENTS serialization system filesystem thread program_options date_time timer chrono regex)

find_package(Boost ${BOOST_FIND_MINIMUM_VERSION} COMPONENTS ${BOOST_FIND_MINIMUM_COMPONENTS})

# Required components
if(NOT Boost_SERIALIZATION_LIBRARY OR NOT Boost_SYSTEM_LIBRARY OR NOT Boost_FILESYSTEM_LIBRARY OR
    NOT Boost_THREAD_LIBRARY OR NOT Boost_DATE_TIME_LIBRARY)
  message(FATAL_ERROR "Missing required Boost components >= v1.65, please install/upgrade Boost or configure your search paths.")
else()
    message("Boost lib found!!")
endif()