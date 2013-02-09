include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(bicycle_printing bicycle_printing.cc)
target_link_libraries(bicycle_printing bicycle)
