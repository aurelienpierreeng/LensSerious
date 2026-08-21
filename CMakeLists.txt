cmake_minimum_required(VERSION 3.18)
project(LensSerious C)

# A library of values: no dependencies beyond libm, on purpose. The parity harness is
# the only thing here that links lensfun, and only when it is found.
add_library(lensserious STATIC src/lensserious.c)
target_include_directories(lensserious PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(lensserious PRIVATE -Wall -Wextra)
if(NOT MSVC)
  target_link_libraries(lensserious PUBLIC m)
endif()
set_target_properties(lensserious PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 11)

option(LENSSERIOUS_TESTS "Build the lensfun parity harness" ON)
if(LENSSERIOUS_TESTS)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(LENSFUN IMPORTED_TARGET lensfun)
  endif()
  if(LENSFUN_FOUND)
    add_executable(parity_lensfun tests/parity_lensfun.c)
    target_link_libraries(parity_lensfun PRIVATE lensserious PkgConfig::LENSFUN)
    # The harness interposes liblensfun's (C++-mangled) CPU detection to force the scalar
    # path; the executable must export its symbols for the interposition to bind.
    target_link_options(parity_lensfun PRIVATE -Wl,--export-dynamic)
    enable_testing()
    add_test(NAME parity COMMAND parity_lensfun)
  else()
    message(STATUS "lensfun not found: the parity harness will not be built here")
  endif()
endif()
