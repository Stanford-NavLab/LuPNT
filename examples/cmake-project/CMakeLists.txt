cmake_minimum_required(VERSION 3.16)

project(app)

find_package(lupnt)

add_executable(app main.cpp)

target_link_libraries(app lupnt::lupnt)
