include_directories(${BICYCLE_SOURCE_DIR}/googletest)
include_directories(${BICYCLE_SOURCE_DIR}/googletest/include)

add_executable(testWhipple testWhipple.cc)
target_link_libraries(testWhipple gtest_main)

# Add test to testing framework
add_test(NAME testWhipple
         COMMAND testWhipple)
