# qcanopy_parser: OpenQASM 3.0 parser (include/qcanopy/parser.hpp).

add_library(qcanopy_parser STATIC
  lexer.cpp
  parser.cpp
)
add_library(qcanopy::parser ALIAS qcanopy_parser)

target_link_libraries(qcanopy_parser PUBLIC qcanopy::ir PRIVATE qcanopy_warnings)

add_executable(qcanopy_parse ${PROJECT_SOURCE_DIR}/tools/parse_files.cpp)
target_link_libraries(qcanopy_parse PRIVATE qcanopy::parser qcanopy_warnings)
