# Unit tests
function(qcanopy_add_test name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE ${ARGN} qcanopy_warnings)
  add_test(NAME ${name} COMMAND ${name})
  # Tests read reference data relative to tests/
  set_tests_properties(${name} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()
 
qcanopy_add_test(test_ir_header cpp/test_ir_header.cpp qcanopy::ir)
 
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cpp/test_gates.cpp)
  qcanopy_add_test(test_gates cpp/test_gates.cpp qcanopy::ir)
endif()