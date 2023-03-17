# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.9)

project(Chrom3D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "" FORCE)
endif()

if(NOT BUILD_SHARED_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(Boost 1.54.0 REQUIRED COMPONENTS filesystem random)

# Try to enable IPO when compiling in Release
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  include(CheckIPOSupported)
  check_ipo_supported(RESULT result OUTPUT output)
  if(result)
    message(STATUS "Enabling Interprocedural optimization")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)

  else()
    message(WARNING "Interprocedural Optimization is not supported. Disabling it!")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
  endif()
endif()

add_executable(Chrom3D)

target_include_directories(Chrom3D PRIVATE src/ src/tclap-1.2.1/include/)

target_sources(
  Chrom3D
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/Bead.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Chrom3D.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Chromosome.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Constraint.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/MCMC.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Model.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Randomizer.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Util.cpp)

target_link_libraries(Chrom3D PUBLIC m Boost::headers Boost::random)

target_compile_features(Chrom3D PRIVATE cxx_std_${CMAKE_CXX_STANDARD})

if(CMAKE_INTERPROCEDURAL_OPTIMIZATION)
  set_target_properties(Chrom3D PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
endif()

set(COMPILER_WARNINGS
    -Wall
    -Wextra
    -Wsign-compare
    -Wunused
    -Wpedantic
    -Wnull-dereference)
target_compile_options(Chrom3D PRIVATE ${COMPILER_WARNINGS})

# Add options to enable/disable sanitizers
option(ENABLE_SANITIZER_ADDRESS "Enable ASAN" OFF)
option(ENABLE_SANITIZER_ADDRESS "Enable LSAN" OFF)
option(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR "Enable UBSAN" OFF)

set(SANITIZERS "")

if(ENABLE_SANITIZER_ADDRESS)
  set(SANITIZERS "address")
endif()

if(ENABLE_SANITIZER_LEAK)
  if("${SANITIZERS}" STREQUAL "")
    set(SANITIZERS "leak")
  else()
    set(SANITIZERS "${SANITIZERS},leak")
  endif()
endif()

if(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR)
  if("${SANITIZERS}" STREQUAL "")
    set(SANITIZERS "undefined")
  else()
    set(SANITIZERS "${SANITIZERS},undefined")
  endif()
endif()

if(NOT
   "${SANITIZERS}"
   STREQUAL
   "")
  target_compile_options(Chrom3D PRIVATE -fsanitize=${SANITIZERS})
  target_link_options(Chrom3D PRIVATE -fsanitize=${SANITIZERS})
endif()

include(GNUInstallDirs)

install(TARGETS Chrom3D BUNDLE DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES LICENCE DESTINATION share/licenses/Chrom3D/)
