
include_directories(${EIGEN3_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_DIRS}
                    ${NUMPY_INCLUDE_DIRS}
                    ${BICYCLE_SOURCE_DIR}/inc)

add_library(bicycle_pywrap 0-example.cpp common.cpp)
target_link_libraries(bicycle_pywrap bicycle ${PYTHON_LIBRARIES})
set_target_properties(bicycle_pywrap PROPERTIES PREFIX "")
