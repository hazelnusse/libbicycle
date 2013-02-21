
include_directories(${EIGEN3_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_DIRS}
                    ${NUMPY_INCLUDE_DIRS})

add_library(bicycle_python 0-example.cpp common.cpp)
target_link_libraries(bicycle_python ${PYTHON_LIBRARIES})
set_target_properties(bicycle_python PROPERTIES PREFIX "")
