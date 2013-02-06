project(BICYCLE CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH
    ${BICYCLE_MODULE_PATH}
    ${BICYCLE_SOURCE_DIR}/cmake/Modules/)

set(BUILD_SHARED_LIBS true)

# Set g++ compiler flags if GNU C++ compiler is being used
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Wall -Wno-unused-local-typedefs -O2 -std=c++11")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wno-unused-local-typedefs -O0 -std=c++11 -ggdb3")
    set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wno-unused-local-typedefs -Weffc++ -O3 -std=c++11")
endif()

find_package(Eigen3 3.1)

include_directories(${BICYCLE_SOURCE_DIR}/inc)
add_subdirectory(${BICYCLE_SOURCE_DIR}/src)
