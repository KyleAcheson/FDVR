cmake_minimum_required(VERSION 3.22)
project(DVR Fortran)

enable_language(Fortran)
set(CMAKE_FC_COMPILER gfortran)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/run)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(CMAKE_PREFIX_PATH /usr/lib/x86_64-linux-gnu/openblas-pthread/cmake/openblas/)
set(OpenBLAS_VERSION "0.3.20")
set(OpenBLAS_INCLUDE_DIRS /usr/include/x86_64-linux-gnu/openblas-pthread/)
set(OpenBLAS_LIBRARIES /usr/lib/x86_64-linux-gnu/openblas-pthread//libopenblas.so)

find_package(OpenBLAS REQUIRED)

file(GLOB SOURCES "src/*.f90")

add_executable(dvr ${SOURCES})

target_compile_options(dvr PRIVATE -O3)
target_link_libraries(dvr ${OpenBLAS_LIBRARIES})
