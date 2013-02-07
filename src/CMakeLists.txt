include_directories(${EIGEN3_INCLUDE_DIR})

add_library(bicycle bicycle.cc bicycle_generated.cc)

if (BICYCLE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

install(TARGETS bicycle
        LIBRARY DESTINATION lib/bicycle
        ARCHIVE DESTINATION lib/bicycle)
