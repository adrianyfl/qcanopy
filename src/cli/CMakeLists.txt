# qcanopy: the user-facing command-line tool (src/cli/main.cpp dispatches to subcommands).
# Subcommand logic is a library so tests can call it without spawning the binary.

add_library(qcanopy_cli STATIC
  check.cpp
)
target_include_directories(qcanopy_cli PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(qcanopy_cli
  PUBLIC qcanopy::parser qcanopy::rewrite qcanopy::equiv
  PRIVATE qcanopy_warnings)
add_library(qcanopy::cli ALIAS qcanopy_cli)

add_executable(qcanopy main.cpp)
target_link_libraries(qcanopy PRIVATE qcanopy::cli qcanopy_warnings)
