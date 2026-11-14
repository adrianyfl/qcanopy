# qcanopy_equiv: simulation-based equivalence checking (include/qcanopy/stimulus.hpp, equiv.hpp).

add_library(qcanopy_equiv STATIC
  stimulus.cpp
)
add_library(qcanopy::equiv ALIAS qcanopy_equiv)

target_link_libraries(qcanopy_equiv PUBLIC qcanopy::sim PRIVATE qcanopy_warnings)
