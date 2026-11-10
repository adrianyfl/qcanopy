# qcanopy_ir: the circuit IR (include/qcanopy/ir.hpp).
 
add_library(qcanopy_ir STATIC
  ir.cpp
  dag.cpp
)
add_library(qcanopy::ir ALIAS qcanopy_ir)

target_include_directories(qcanopy_ir PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_compile_features(qcanopy_ir PUBLIC cxx_std_17)
target_link_libraries(qcanopy_ir PRIVATE qcanopy_warnings)
 