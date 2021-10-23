cmake_minimum_required(VERSION 3.20)
project(MonteCarlo_option_obj)

set(CMAKE_CXX_STANDARD 14)

add_executable(MonteCarlo_option_obj main.cpp utility.h utility.cpp rv_library.h rv_library.cpp EuroCall.h EuroCall.cpp GBM_process.h GBM_process.cpp Accumulator.cpp Accumulator.h)
