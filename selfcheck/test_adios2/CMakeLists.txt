cmake_minimum_required(VERSION 3.12)
project(mini C CXX)

find_package(ADIOS2 REQUIRED)

add_executable(mini mini.cpp)
target_link_libraries(mini PRIVATE adios2::cxx11)