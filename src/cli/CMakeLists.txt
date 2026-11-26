# qcanopy: the user-facing command-line tool. main.cpp dispatches to the subcommands in
# commands.hpp: check (the whole pipeline), parse, rewrite, sim, generate.
# Subcommand logic is a library so tests, and the older qcanopy_parse / qcanopy_generate
# executables (built in src/parser and src/generator), can use it without the binary.
#
# `check` runs the Python pair prover in rewriting/ as a child process (rewriter.hpp). The
# interpreter and run.py directory found here are defaults; QCANOPY_PYTHON and
# QCANOPY_REWRITER override them at run time.

add_library(qcanopy_cli STATIC
  common.cpp
  check.cpp
  parse.cpp
  rewrite.cpp
  rewriter.cpp
  sim.cpp
  generate.cpp
)
target_include_directories(qcanopy_cli
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
  PRIVATE ${PROJECT_SOURCE_DIR}/src/generator)   # json.hpp (private strict JSON reader)
target_link_libraries(qcanopy_cli
  PUBLIC qcanopy::parser qcanopy::equiv
  PRIVATE qcanopy::generator qcanopy::rewrite qcanopy_warnings)   # dataset I/O; normalize, peel

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
target_compile_definitions(qcanopy PRIVATE QCANOPY_VERSION="${PROJECT_VERSION}")
