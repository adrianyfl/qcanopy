# Unit tests
function(qcanopy_add_test name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE ${ARGN} qcanopy_warnings)
  add_test(NAME ${name} COMMAND ${name})
  # Tests read reference data relative to tests/
  set_tests_properties(${name} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()
 
qcanopy_add_test(test_ir_header cpp/test_ir_header.cpp qcanopy::ir)
qcanopy_add_test(test_lexer cpp/test_lexer.cpp qcanopy::parser)
qcanopy_add_test(test_parser cpp/test_parser.cpp qcanopy::parser)
qcanopy_add_test(test_parser_roundtrip cpp/test_parser_roundtrip.cpp qcanopy::parser qcanopy::generator)
qcanopy_add_test(test_normalize cpp/test_normalize.cpp qcanopy::rewrite)
qcanopy_add_test(test_generator cpp/test_generator.cpp qcanopy::generator)
qcanopy_add_test(test_dataset_io cpp/test_dataset_io.cpp qcanopy::generator)
qcanopy_add_test(test_qasm_writer cpp/test_qasm_writer.cpp qcanopy::generator)
qcanopy_add_test(test_stimulus cpp/test_stimulus.cpp qcanopy::equiv)
qcanopy_add_test(test_equiv cpp/test_equiv.cpp qcanopy::equiv)
qcanopy_add_test(test_cli_check cpp/test_cli_check.cpp qcanopy::cli)
if(Python3_Interpreter_FOUND)
  target_compile_definitions(test_cli_check PRIVATE QCANOPY_TEST_HAVE_PYTHON=1)
else()
  target_compile_definitions(test_cli_check PRIVATE QCANOPY_TEST_HAVE_PYTHON=0)
endif()
add_test(NAME test_check_cli COMMAND ${CMAKE_COMMAND}
  "-DQCANOPY=$<TARGET_FILE:qcanopy>" "-DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}"
  -P "${CMAKE_CURRENT_SOURCE_DIR}/test_check_cli.cmake")
if(Python3_Interpreter_FOUND)
  set(QCANOPY_TEST_HAVE_PYTHON ON)
else()
  set(QCANOPY_TEST_HAVE_PYTHON OFF)
endif()
add_test(NAME test_cli_commands COMMAND ${CMAKE_COMMAND}
  "-DQCANOPY=$<TARGET_FILE:qcanopy>" "-DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}"
  "-DHAVE_PYTHON=${QCANOPY_TEST_HAVE_PYTHON}"
  -P "${CMAKE_CURRENT_SOURCE_DIR}/test_cli_commands.cmake")
add_test(NAME test_generator_cli COMMAND ${CMAKE_COMMAND}
  "-DGENERATOR=$<TARGET_FILE:qcanopy_generate>" "-DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}"
  -P "${CMAKE_CURRENT_SOURCE_DIR}/test_cli.cmake")
 
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cpp/test_gates.cpp)
  qcanopy_add_test(test_gates cpp/test_gates.cpp qcanopy::ir)
endif()
