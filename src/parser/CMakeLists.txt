# qcanopy_parser: circuit source parser (include/qcanopy/parser.hpp).

add_library(qcanopy_parser STATIC
  lexer.cpp
  parser.cpp
)
add_library(qcanopy::parser ALIAS qcanopy_parser)

target_link_libraries(qcanopy_parser PUBLIC qcanopy::ir PRIVATE qcanopy_warnings)
