include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(constraint_forces constraint_forces.cc)
target_link_libraries(constraint_forces bicycle)

add_executable(dependent_speeds dependent_speeds.cc)
target_link_libraries(dependent_speeds bicycle)

add_executable(stream_output stream_output.cc)
target_link_libraries(stream_output bicycle)
