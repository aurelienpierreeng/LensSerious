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

option(LENSSERIOUS_TESTS "Build the lensfun parity harness" ON)
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
    add_test(NAME parity COMMAND parity_lensfun)

    # The importer reads the XML through liblensfun on purpose -- see the file header. It
    # is an offline tool, so that dependency never reaches a render.
    if(TARGET lensserious_db)
      add_executable(import_lensfun_xml tools/import_lensfun_xml.c)
      target_link_libraries(import_lensfun_xml PRIVATE lensserious_db PkgConfig::LENSFUN)
      target_compile_options(import_lensfun_xml PRIVATE -Wall -Wextra)
      set_target_properties(import_lensfun_xml PROPERTIES C_STANDARD 11)

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
      set_tests_properties(db_match PROPERTIES FIXTURES_REQUIRED LENSDB)
      set_tests_properties(db_import PROPERTIES FIXTURES_SETUP LENSDB)
      set_tests_properties(db_parity PROPERTIES FIXTURES_REQUIRED LENSDB)
    endif()
  else()
    message(STATUS "lensfun not found: the parity harness will not be built here")
  endif()
endif()
