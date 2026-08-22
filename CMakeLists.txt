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

# The database is an OPTIONAL component: the evaluators above stay dependency-free, since
# they are what ends up inside a pixel pipeline. Only the lookup half needs SQLite.
option(LENSSERIOUS_DB "Build the SQLite database reader and importer" ON)
if(LENSSERIOUS_DB)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(SQLITE3 IMPORTED_TARGET sqlite3)
  endif()
  if(SQLITE3_FOUND)
    add_library(lensserious_db STATIC src/lensserious_db.c)
    target_include_directories(lensserious_db PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(lensserious_db PUBLIC lensserious PkgConfig::SQLITE3)
    target_compile_options(lensserious_db PRIVATE -Wall -Wextra)
    set_target_properties(lensserious_db PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 11)
  else()
    message(STATUS "sqlite3 not found: the database reader will not be built")
  endif()
endif()

# liblensfun, found once and used by two different things below: the offline importer,
# which is a TOOL a consumer needs in order to have a database at all, and the parity
# harnesses, which are tests. The importer used to live inside the test block, which meant
# -DLENSSERIOUS_TESTS=OFF silently removed the only way to build a database -- exactly what
# a consumer embedding this library wants to turn off.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(LENSFUN IMPORTED_TARGET lensfun)
endif()

# The importer reads the XML through liblensfun on purpose -- see the file header. It is an
# offline tool run when upstream publishes calibrations, so that dependency never reaches a
# render: nothing the pixel pipeline links to knows lensfun exists.
if(LENSFUN_FOUND AND TARGET lensserious_db)
  # lf_db_load_directory() arrived in liblensfun 0.3.3. Ubuntu 22.04 -- what this project's
  # consumers build their nightlies on -- still ships 0.3.2, where the symbol is absent and
  # the link fails with an undefined reference rather than anything that names a version.
  # Probed rather than inferred from the version string, because a distribution's version is
  # not a reliable statement about which symbols it exports.
  include(CheckSymbolExists)
  set(CMAKE_REQUIRED_INCLUDES ${LENSFUN_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES PkgConfig::LENSFUN)
  check_symbol_exists(lf_db_load_directory "lensfun.h" LS_HAVE_LF_DB_LOAD_DIRECTORY)
  unset(CMAKE_REQUIRED_INCLUDES)
  unset(CMAKE_REQUIRED_LIBRARIES)

  add_executable(import_lensfun_xml tools/import_lensfun_xml.c)
  target_link_libraries(import_lensfun_xml PRIVATE lensserious_db PkgConfig::LENSFUN)
  target_compile_options(import_lensfun_xml PRIVATE -Wall -Wextra)
  set_target_properties(import_lensfun_xml PROPERTIES C_STANDARD 11)
  if(LS_HAVE_LF_DB_LOAD_DIRECTORY)
    target_compile_definitions(import_lensfun_xml PRIVATE LS_HAVE_LF_DB_LOAD_DIRECTORY)
  else()
    message(STATUS "liblensfun has no lf_db_load_directory (pre-0.3.3): "
                   "the importer will walk the directory itself")
  endif()
endif()

option(LENSSERIOUS_TESTS "Build the lensfun parity harness" ON)
# Every parity, oracle and benchmark target needs liblensfun, and the whole point of this
# repository's CI is that they RUN. Finding lensfun is therefore optional for a developer
# (who may not have it) and mandatory in CI: without this, a missing dependency degrades to
# "configured fine, built nothing", and the failure surfaces later as a missing binary
# rather than as the missing package it actually is.
option(LENSSERIOUS_REQUIRE_LENSFUN "Fail configuration if liblensfun is not found" OFF)
if(LENSSERIOUS_TESTS)
  enable_testing()

  # include/lensserious_eval.h is compiled twice: as C99 into the library above, and as
  # OpenCL C by whatever driver the host has at run time. Only the first of those is
  # covered by building this project, so the second rots quietly -- which is how the old
  # hand-copied kernel came to use native_sqrt() and to have no vignetting evaluator at
  # all. clang compiles OpenCL C offline, without a device, so the second compilation can
  # be a test like any other.
  find_program(LENSSERIOUS_CLANG NAMES clang)
  if(LENSSERIOUS_CLANG)
    add_test(NAME opencl_syntax
             COMMAND ${LENSSERIOUS_CLANG}
                     -x cl -cl-std=CL1.2 -Xclang -finclude-default-header
                     --target=spir -fsyntax-only -Werror
                     -I ${CMAKE_CURRENT_SOURCE_DIR}/include
                     ${CMAKE_CURRENT_SOURCE_DIR}/opencl/lensserious.cl)
  else()
    message(STATUS "clang not found: the OpenCL kernel will not be syntax-checked here")
  endif()

  # The maker-profile path is the one part of this library liblensfun cannot arbitrate --
  # it has no knot model to compare against -- so this test asserts properties instead, and
  # is deliberately outside the LENSFUN_FOUND block: it must run everywhere, including on
  # the platforms where liblensfun is not installed at all.
  add_executable(knots tests/knots.c)
  target_link_libraries(knots PRIVATE lensserious)
  target_compile_options(knots PRIVATE -Wall -Wextra)
  set_target_properties(knots PROPERTIES C_STANDARD 11)
  add_test(NAME knots COMMAND knots)

  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(LENSFUN IMPORTED_TARGET lensfun)
  endif()
  if(NOT LENSFUN_FOUND AND LENSSERIOUS_REQUIRE_LENSFUN)
    message(FATAL_ERROR
            "liblensfun was requested (LENSSERIOUS_REQUIRE_LENSFUN=ON) but pkg-config "
            "cannot find it. Note that lensfun.pc requires glib-2.0, so libglib2.0-dev "
            "must be installed as well -- without it the module resolves as not-found "
            "with no mention of lensfun itself.")
  endif()
  if(LENSFUN_FOUND)
    add_executable(parity_lensfun tests/parity_lensfun.c)
    target_link_libraries(parity_lensfun PRIVATE lensserious PkgConfig::LENSFUN)
    # The harness interposes liblensfun's (C++-mangled) CPU detection to force the scalar
    # path; the executable must export its symbols for the interposition to bind.
    target_link_options(parity_lensfun PRIVATE -Wl,--export-dynamic)
    add_test(NAME parity COMMAND parity_lensfun)

    if(TARGET import_lensfun_xml)
      add_executable(parity_db tests/parity_db.c)
      target_link_libraries(parity_db PRIVATE lensserious_db PkgConfig::LENSFUN)
      target_compile_options(parity_db PRIVATE -Wall -Wextra)
      set_target_properties(parity_db PROPERTIES C_STANDARD 11)
      # Rebuild the database from the system XML, then hold it to liblensfun's own answer.
      add_test(NAME db_import
               COMMAND import_lensfun_xml ${CMAKE_CURRENT_SOURCE_DIR}/src/schema.sql
                       ${CMAKE_CURRENT_BINARY_DIR}/lenses.db)
      add_test(NAME db_parity COMMAND parity_db ${CMAKE_CURRENT_BINARY_DIR}/lenses.db)

      add_executable(match_lensfun tests/match_lensfun.c)
      target_link_libraries(match_lensfun PRIVATE lensserious_db PkgConfig::LENSFUN)
      target_compile_options(match_lensfun PRIVATE -Wall -Wextra)
      set_target_properties(match_lensfun PROPERTIES C_STANDARD 11)
      add_test(NAME db_match COMMAND match_lensfun ${CMAKE_CURRENT_BINARY_DIR}/lenses.db)

      # The OpenCL half of the comparison: what a correction costs to REACH a GPU, which
      # is where this library's claim actually lives. Optional -- it needs a device.
      find_path(OPENCL_INCLUDE_DIR CL/cl.h)
      # libOpenCL.so.1 without a devel symlink is the normal state on a machine that has a
      # driver but no ICD-loader -devel package, so accept the versioned name too.
      find_library(OPENCL_LIBRARY NAMES OpenCL libOpenCL.so.1)
      if(OPENCL_INCLUDE_DIR AND OPENCL_LIBRARY)
        add_executable(bench_opencl tests/bench_opencl.c)
        target_include_directories(bench_opencl PRIVATE ${OPENCL_INCLUDE_DIR})
        target_link_libraries(bench_opencl PRIVATE lensserious_db PkgConfig::LENSFUN
                              ${OPENCL_LIBRARY} m)
        target_compile_options(bench_opencl PRIVATE -Wall -Wextra)
        target_compile_definitions(bench_opencl PRIVATE
          LS_EVAL_HEADER_PATH="${CMAKE_CURRENT_SOURCE_DIR}/include/lensserious_eval.h"
          LS_KERNEL_PATH="${CMAKE_CURRENT_SOURCE_DIR}/opencl/lensserious.cl")
        set_target_properties(bench_opencl PROPERTIES C_STANDARD 11)
      else()
        message(STATUS "OpenCL not found: the device benchmark will not be built")
      endif()

      # Not a test: it measures, it does not assert. Run it by hand.
      add_executable(bench_lensfun tests/bench_lensfun.c)
      target_link_libraries(bench_lensfun PRIVATE lensserious_db PkgConfig::LENSFUN m)
      # A pixel pipeline runs these stages across threads; measure them that way too.
      find_package(OpenMP COMPONENTS C)
      if(OpenMP_C_FOUND)
        target_link_libraries(bench_lensfun PRIVATE OpenMP::OpenMP_C)
      endif()
      target_compile_options(bench_lensfun PRIVATE -Wall -Wextra)
      set_target_properties(bench_lensfun PROPERTIES C_STANDARD 11)
      set_tests_properties(db_match PROPERTIES FIXTURES_REQUIRED LENSDB)
      set_tests_properties(db_import PROPERTIES FIXTURES_SETUP LENSDB)
      set_tests_properties(db_parity PROPERTIES FIXTURES_REQUIRED LENSDB)
    endif()
  else()
    message(STATUS "lensfun not found: the parity harness will not be built here")
  endif()
endif()
