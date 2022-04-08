
cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(dune_pisces VERSION 00.01 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
#find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

# Find our own CMake modules, and tell everyone else where to find them.
cet_cmake_module_directories(Modules BINARY)

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-expansion-to-defined
)

cet_report_compiler_flags()

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(Core)

# ups - table and config files
add_subdirectory(ups)

# Dependencies
include_directories(
  $ENV{PISCESOSC_INC}
)
