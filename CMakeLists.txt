cmake_minimum_required (VERSION 3.13.0)
project (GeographicLib-Documentation NONE)

find_program (SPHINX sphinx-build REQUIRED)
find_program (RSYNC rsync)

add_subdirectory (geodesic-papers)
add_subdirectory (misc)
add_subdirectory (source)

set (BINARY_EXT "gif|pdf|png|kmz")
add_custom_target (checktrailingspace
  COMMAND git ls-files |
  grep -E -v '\\.\(${BINARY_EXT}\)$$' |
  xargs grep '[ \t]$$' || true
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMENT "Looking for trailing spaces")
add_custom_target (checktabs
  COMMAND git ls-files |
  grep -E -v '\\.\(${BINARY_EXT}\)$$' |
  xargs grep '\t' || true
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMENT "Looking for tabs")
add_custom_target (checkblanklines
  COMMAND git ls-files |
  grep -E -v '\\.\(${BINARY_EXT}\)$$' |
  while read f\; do tr 'X\\n' 'YX' < $$f |
  grep -E '\(^X|XXX|XX$$|[^X]$$\)' > /dev/null && echo $$f\; done || true
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMENT "Looking for extra blank lines")

add_custom_target (sanitize)
add_dependencies (sanitize checktrailingspace checktabs checkblanklines)
