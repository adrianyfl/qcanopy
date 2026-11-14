# Unit tests
function(qcanopy_add_test name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE ${ARGN} qcanopy_warnings)
  add_test(NAME ${name} COMMAND ${name})
  # Tests read reference data relative to tests/
  set_tests_properties(${name} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()
 
qcanopy_add_test(test_ir_header cpp/test_ir_header.cpp qcanopy::ir)
qcanopy_add_test(test_normalize cpp/test_normalize.cpp qcanopy::rewrite)
qcanopy_add_test(test_generator cpp/test_generator.cpp qcanopy::generator)
qcanopy_add_test(test_dataset_io cpp/test_dataset_io.cpp qcanopy::generator)
qcanopy_add_test(test_qasm_writer cpp/test_qasm_writer.cpp qcanopy::generator)
qcanopy_add_test(test_stimulus cpp/test_stimulus.cpp qcanopy::equiv)
add_test(NAME test_generator_cli COMMAND ${CMAKE_COMMAND}
  "-DGENERATOR=$<TARGET_FILE:qcanopy_generate>" "-DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}"
  -P "${CMAKE_CURRENT_SOURCE_DIR}/test_cli.cmake")
 
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cpp/test_gates.cpp)
  qcanopy_add_test(test_gates cpp/test_gates.cpp qcanopy::ir)
endif()
