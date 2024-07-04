# The path where cmake config files are installed
set(LUPNT_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/lupnt)

install(EXPORT lupntTargets
    FILE lupntTargets.cmake
    NAMESPACE lupnt::
    DESTINATION ${LUPNT_INSTALL_CONFIGDIR}
    COMPONENT cmake)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/lupntConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    ARCH_INDEPENDENT)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/lupntConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/lupntConfig.cmake
    INSTALL_DESTINATION ${LUPNT_INSTALL_CONFIGDIR}
    PATH_VARS LUPNT_INSTALL_CONFIGDIR)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/lupntConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/lupntConfigVersion.cmake
    DESTINATION ${LUPNT_INSTALL_CONFIGDIR})
