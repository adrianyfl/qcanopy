# qcanopy_rewrite: normalization and rewrite passes over CircuitDAG (include/qcanopy/normalize.hpp).

add_library(qcanopy_rewrite STATIC
  normalize.cpp
  peel.cpp
)
add_library(qcanopy::rewrite ALIAS qcanopy_rewrite)

target_link_libraries(qcanopy_rewrite PUBLIC qcanopy::ir PRIVATE qcanopy_warnings)
