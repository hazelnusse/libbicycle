include_directories(${BICYCLE_SOURCE_DIR}/googletest)
include_directories(${BICYCLE_SOURCE_DIR}/googletest/include)

add_executable(testWhipple testWhipple.cc)
target_link_libraries(testWhipple gtest_main)

add_executable(testNormalForces testNormalForces.cc)
target_link_libraries(testNormalForces gtest_main bicycle)

add_executable(testKinematics testKinematics.cc)
target_link_libraries(testKinematics gtest_main bicycle)

# Add test to testing framework
add_test(NAME testWhipple
         COMMAND testWhipple)

add_test(NAME testNormalForces
         COMMAND testNormalForces)

add_test(NAME testKinematics
         COMMAND testKinematics)
