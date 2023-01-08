cmake_minimum_required(VERSION 3.23)
project(mda)

set(CMAKE_CXX_STANDARD 14)

add_executable(mda main.cpp Equation.cpp Equation.h Grid.cpp Grid.h Solver.cpp Solver.h FDM.cpp FDM.h GalerkinMethod.cpp GalerkinMethod.h const.h)
