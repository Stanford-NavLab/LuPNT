file(GLOB_RECURSE CPPFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cc)

find_package(OpenMP REQUIRED)

foreach(CPPFILE ${CPPFILES})
    get_filename_component(CPPNAME ${CPPFILE} NAME_WE)
    get_filename_component(CPPDIR ${CPPFILE} DIRECTORY)

    add_executable(${CPPNAME} ${CPPFILE})
    target_link_libraries(${CPPNAME} PRIVATE lupnt::lupnt)
    set_target_properties(${CPPNAME} PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
    target_link_libraries(${CPPNAME} PRIVATE OpenMP::OpenMP_CXX)
    # Set target folder for IDEs that support it, under "examples"
    set_target_properties(${CPPNAME} PROPERTIES FOLDER "examples/${CPPDIR}")
endforeach()