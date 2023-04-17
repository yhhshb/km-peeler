cmake_minimum_required(VERSION 3.11)
project(KMP)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (UNIX AND NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2") # for hardware popcount and pdep
endif()

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-char-subscripts")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  if (KMP_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()
endif()

include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

set(KMP_SOURCES
  bundled/biolib/src/io.cpp
  bundled/biolib/src/hash.cpp
  bundled/biolib/bundled/MurmurHash3.cpp
  src/build.cpp
  src/diff.cpp
  src/correct.cpp
  src/dump.cpp
  src/bias.cpp
  src/IBLT.cpp
  src/utils.cpp
)

set(CWS_SOURCES
  bundled/biolib/bundled/MurmurHash3.cpp
  scripts/kmc_extra/kmc_api/kmc_file.cpp
  scripts/kmc_extra/kmc_api/kmer_api.cpp
  scripts/kmc_extra/kmc_api/mmer.cpp
)

add_executable(kmp src/kmp.cpp ${KMP_SOURCES})
target_link_libraries(kmp argparse z)
add_executable(cws scripts/kmc_extra/cws.cpp ${CWS_SOURCES})
target_link_libraries(cws argparse)