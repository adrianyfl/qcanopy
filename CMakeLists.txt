cmake_minimum_required(VERSION 3.18)
project(qcanopy VERSION 0.1.0 LANGUAGES CXX)
 
option(QCANOPY_BUILD_TESTS "Build unit tests" ON)
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 
 
# Shared warning flags: link into every target we own
add_library(qcanopy_warnings INTERFACE)
if(MSVC)
  target_compile_options(qcanopy_warnings INTERFACE /W4 /permissive-)
else()
  target_compile_options(qcanopy_warnings INTERFACE -Wall -Wextra -Wpedantic)
endif()
 
# Modules
add_subdirectory(src/ir)                 # qcanopy::ir
add_subdirectory(src/rewrite)            # qcanopy::rewrite
add_subdirectory(src/sim)                # qcanpoy::sim
add_subdirectory(src/generator)          # qcanopy::generator, qcanopy_generate
 
# Tests
if(QCANOPY_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
