include_directories(${EIGEN3_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_DIRS}
                    ${NUMPY_INCLUDE_DIRS}
                    ${BICYCLE_SOURCE_DIR}/inc
                    ${Boost_INCLUDE_DIRS})

add_library(whipple_ext whipple_ext.cc)
target_link_libraries(whipple_ext ${Boost_LIBRARIES})
set_target_properties(whipple_ext PROPERTIES
                      PREFIX ""
                      INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

add_library(wheelassemblygyrostat_ext wheelassemblygyrostat_ext.cc)
target_link_libraries(wheelassemblygyrostat_ext bicycle ${Boost_LIBRARIES})
set_target_properties(wheelassemblygyrostat_ext PROPERTIES
                      PREFIX ""
                      INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

add_library(bicycle_ext bicycle_ext.cc)
target_link_libraries(bicycle_ext bicycle ${Boost_LIBRARIES})
set_target_properties(bicycle_ext PROPERTIES
                      PREFIX ""
                      INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

install(TARGETS bicycle_ext wheelassemblygyrostat_ext whipple_ext
        DESTINATION lib64/python3.2/site-packages/bicycle)
install(FILES __init__.py
        DESTINATION lib64/python3.2/site-packages/bicycle)
