# qcanopy: the user-facing command-line tool (src/cli/main.cpp dispatches to subcommands).
# Subcommand logic is a library so tests can call it without spawning the binary.
#
# `check` runs the Python pair prover in rewriting/ as a child process (rewriter.hpp). The
# interpreter and run.py directory found here are defaults; QCANOPY_PYTHON and
# QCANOPY_REWRITER override them at run time.

add_library(qcanopy_cli STATIC
  check.cpp
  rewriter.cpp
)
target_include_directories(qcanopy_cli
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
  PRIVATE ${PROJECT_SOURCE_DIR}/src/generator)   # json.hpp (private strict JSON reader)
target_link_libraries(qcanopy_cli
  PUBLIC qcanopy::parser qcanopy::equiv
  PRIVATE qcanopy::generator qcanopy_warnings)   # write_ir

if(Python3_Interpreter_FOUND)
  set(QCANOPY_PYTHON_DEFAULT "${Python3_EXECUTABLE}")
else()
  set(QCANOPY_PYTHON_DEFAULT "python3")
endif()
target_compile_definitions(qcanopy_cli PRIVATE
  QCANOPY_REWRITER_DIR="${PROJECT_SOURCE_DIR}/rewriting"
  QCANOPY_PYTHON_DEFAULT="${QCANOPY_PYTHON_DEFAULT}")
add_library(qcanopy::cli ALIAS qcanopy_cli)

add_executable(qcanopy main.cpp)
target_link_libraries(qcanopy PRIVATE qcanopy::cli qcanopy_warnings)
