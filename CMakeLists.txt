project (GeographicLib-Documentation NONE)

cmake_minimum_required (VERSION 3.13.0)

find_program (SPHINX sphinx-build REQUIRED)
find_program (RSYNC rsync)

add_subdirectory (geodesic-papers)
add_subdirectory (misc)
add_subdirectory (source)
