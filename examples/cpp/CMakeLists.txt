#### Create Examples

add_executable(example1 ${CMAKE_CURRENT_SOURCE_DIR}/example1.cpp)
target_link_libraries(example1 PRIVATE LuPNT::CPPL)