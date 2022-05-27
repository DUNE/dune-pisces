
cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)

project(DUNE_PISCES VERSION 1.00)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

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

# dependencies
#find_ups_product(cafanacore)
find_ups_product(dunecore)

list(APPEND CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -DTBB_INTERFACE_NEW)

include_directories(
  $ENV{PISCESOSC_INC}
  $ENV{ROOT_INC}
  $ENV{CAFANA_INC}
  $ENV{CAFANACORE_INC}
  $ENV{STAN_MATH_INC}
  $ENV{STAN_INC}
  $ENV{EIGEN_INC}
  $ENV{TBB_INC}
  $ENV{DUNEPISCES_DIR}
)

function(link_root target)
    unset(rootlibs CACHE)
    list(APPEND rootlibs ${ROOT_LIBRARIES})
    foreach(lib ${rootlibs})
        target_link_libraries(${target} PUBLIC ${lib})
    endforeach()
endfunction()

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(Core)

# ups - table and config files
add_subdirectory(ups)
