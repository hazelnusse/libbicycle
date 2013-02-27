include_directories(${EIGEN3_INCLUDE_DIR})

add_library(bicycle bicycle_convert_whipple.cc
                    bicycle_generated.cc
                    bicycle_interface.cc
                    bicycle_kinematics.cc
                    bicycle_linearization.cc
                    bicycle_numerical.cc
                    wheelassemblygyrostat.cc)

if (BICYCLE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

install(TARGETS bicycle
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
