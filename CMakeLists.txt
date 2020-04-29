cmake_minimum_required (VERSION 3.7.2)

project (nemss2)

set (nemss2_VERSION_MAJOR 0)
set (nemss2_VERSION_MAJOR 1)

# googletest download & unpack

configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
  ${CMAKE_BINARY_DIR}/googletest-build
  EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

# Your build code here

include(CTest)
enable_testing()

find_library(LIB_BPP_CORE bpp-core HINTS "/home/aritchie/local/bpp/dev" NO_SYSTEM_ENVIONMENT_PATH)
find_library(LIB_BPP_PHYL bpp-phyl HINTS "/home/aritchie/local/bpp/dev" NO_SYSTEM_ENVIRONMENT_PATH)
find_library(LIB_BPP_SEQ bpp-seq HINTS "/home/aritchie/local/bpp/dev" NO_SYSTEM_ENVIRONMENT_PATH)

set(BPP_LIBS ${LIB_BPP_CORE} ${LIB_BPP_SEQ} ${LIB_BPP_PHYL})
set(CMAKE_CXX_STANDARD 11)

configure_file(nemss_config.h.in ${CMAKE_BINARY_DIR}/nemss_config.h)
include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(src)

set(BUILD_TESTING False)

if(BUILD_TESTING)
  add_subdirectory(test)
  add_test(NAME basictest COMMAND basictest WORKING_DIRECTORY test)
  add_test(NAME mutseltest COMMAND mutseltest WORKING_DIRECTORY test)
  set_tests_properties(basictest mutseltest PROPERTIES ENVIRONMENT LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/aritchie/local/bpp/dev)
endif()
