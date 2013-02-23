include_directories(${EIGEN3_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_DIRS}
                    ${NUMPY_INCLUDE_DIRS}
                    ${BICYCLE_SOURCE_DIR}/inc
                    ${Boost_INCLUDE_DIRS})

add_library(whipple_ext whipple_ext.cc)
target_link_libraries(whipple_ext ${Boost_LIBRARIES})
set_target_properties(whipple_ext PROPERTIES PREFIX "")

add_library(wheelassemblygyrostat_ext wheelassemblygyrostat_ext.cc)
target_link_libraries(wheelassemblygyrostat_ext bicycle ${Boost_LIBRARIES})
set_target_properties(wheelassemblygyrostat_ext PROPERTIES PREFIX "")

add_library(bicycle_ext bicycle_ext.cc)
target_link_libraries(bicycle_ext bicycle ${Boost_LIBRARIES})
set_target_properties(bicycle_ext PROPERTIES PREFIX "")
