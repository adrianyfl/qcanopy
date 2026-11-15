# qcanopy_equiv: simulation-based equivalence checking (include/qcanopy/equiv.hpp, stimulus.hpp).
# The checker owns a qcanopy::sim Simulator and drives it with stimulus batches.

add_library(qcanopy_equiv STATIC
  equiv.cpp
  stimulus.cpp
)
add_library(qcanopy::equiv ALIAS qcanopy_equiv)

target_link_libraries(qcanopy_equiv PUBLIC qcanopy::sim PRIVATE qcanopy_warnings)
