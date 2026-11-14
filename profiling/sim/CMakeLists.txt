# sim_bench: worker executable for profiling/sim/sim_profile.py.
# Add to the top-level CMakeLists.txt with:  add_subdirectory(profiling/sim)

add_executable(sim_bench sim_bench.cpp)
target_compile_features(sim_bench PRIVATE cxx_std_17)
target_include_directories(sim_bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# The simulator: use an existing library target if the sim branch defines one,
# otherwise compile sim.cpp directly into the bench.
set(QCANOPY_SIM_SOURCE "${PROJECT_SOURCE_DIR}/src/sim/sim.cpp" CACHE FILEPATH
    "sim.cpp to compile when no simulator library target exists")
if(TARGET qcanopy::sim)
  target_link_libraries(sim_bench PRIVATE qcanopy::sim)
elseif(TARGET qcanopy_sim)
  target_link_libraries(sim_bench PRIVATE qcanopy_sim)
elseif(EXISTS "${QCANOPY_SIM_SOURCE}")
  target_sources(sim_bench PRIVATE "${QCANOPY_SIM_SOURCE}")
else()
  message(FATAL_ERROR "sim_bench: no qcanopy::sim target and no ${QCANOPY_SIM_SOURCE}; "
                      "set QCANOPY_SIM_SOURCE to the path of sim.cpp")
endif()

# IR, DAG, read_ir and normalized_copy all come through the generator library.
target_link_libraries(sim_bench PRIVATE qcanopy::generator)

if(WIN32)
  target_link_libraries(sim_bench PRIVATE psapi)
endif()

# Profiling numbers are only meaningful from an optimized build.
if(CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE MATCHES "Release|RelWithDebInfo")
  message(WARNING "sim_bench: CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}; use Release for timing")
endif()
