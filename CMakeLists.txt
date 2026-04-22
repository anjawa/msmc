cmake_minimum_required(VERSION 3.15)
project(msmc_cpp LANGUAGES CXX)

find_package(
  Python 3.11 REQUIRED
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(
  _core
  STABLE_ABI
  NB_STATIC
  src/cpp/bindings.cpp
  src/cpp/mayer/pair_potential.cpp
  src/cpp/mayer/hard_sphere.cpp
  src/cpp/mayer/lennard_jones.cpp
  src/cpp/mcmc/metropolis_hasting.cpp)
target_include_directories(_core PRIVATE src/cpp)
install(TARGETS _core LIBRARY DESTINATION msmc)
