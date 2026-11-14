add_library(qcanopy_sim STATIC sim.cpp)
add_library(qcanopy::sim ALIAS qcanopy_sim)

target_include_directories(qcanopy_sim PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_compile_features(qcanopy_sim PUBLIC cxx_std_17)
target_link_libraries(qcanopy_sim PUBLIC qcanopy::ir)
target_link_libraries(qcanopy_sim PRIVATE qcanopy_warnings)