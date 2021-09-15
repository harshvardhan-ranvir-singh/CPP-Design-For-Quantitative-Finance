cmake_minimum_required(VERSION 3.19)
project(MC_v0_plain)

set(CMAKE_CXX_STANDARD 11)

add_executable(MC_v0_plain main.cpp utility.h utility.cpp rv_library.h rv_library.cpp)