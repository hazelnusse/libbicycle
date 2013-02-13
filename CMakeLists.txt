project(BICYCLE CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${BICYCLE_SOURCE_DIR}/cmake/Modules)

option(BUILD_SHARED_LIBS "Build shared libraries." ON)
option(BICYCLE_BUILD_TESTS "Build tests." ON)
option(BICYCLE_BUILD_EXAMPLES "Build examples." ON)
option(BICYCLE_BUILD_DOXYGEN_DOCS "Build doxygen documentation." ON)
option(BICYCLE_BUILD_PYTHON_EXTENSION_MODULES "Build extension modules for Python" ON)

# Set g++ compiler flags if GNU C++ compiler is being used
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Wall -Wno-unused-local-typedefs -O2 -std=c++11")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wno-unused-local-typedefs -O0 -std=c++11 -ggdb3")
    set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wno-unused-local-typedefs -O3 -std=c++11")
endif()

find_package(Eigen3 3.1)

if (BICYCLE_BUILD_TESTS)
  enable_testing()
  add_subdirectory(googletest)
endif()

include_directories(${BICYCLE_SOURCE_DIR}/inc)
add_subdirectory(${BICYCLE_SOURCE_DIR}/inc)
add_subdirectory(${BICYCLE_SOURCE_DIR}/src)

if (BICYCLE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# If we find git installed, set OBD_VERSION_COMMIT to most recent SHA1
find_package(Git)
if(GIT_FOUND)
  # TODO:  check that ${PROJECT_SOURCE_DIR} has a .git folder to ensure it is
  # a git repo rather than just a source tar-ball.
  execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --pretty=oneline
                  WORKING_DIRECTORY ${BICYCLE_SOURCE_DIR}
                  OUTPUT_VARIABLE BICYCLE_VERSION_COMMIT)

  # Keep everything up to the first space, add quotes
  string(REGEX REPLACE "([^ ]+).*"
         "\"\\1\"" BICYCLE_VERSION_COMMIT
         ${BICYCLE_VERSION_COMMIT})
else()
  set(BICYCLE_VERSION_COMMIT "\"\"")
endif()

if (BICYCLE_BUILD_DOXYGEN_DOCS)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
      # prepare doxygen configuration file
      configure_file(${BICYCLE_SOURCE_DIR}/Doxyfile.in
                     ${BICYCLE_BINARY_DIR}/Doxyfile)
      add_custom_target(doc
                        ${DOXYGEN_EXECUTABLE}
                        ${BICYCLE_BINARY_DIR}/Doxyfile
                        WORKING_DIRECTORY ${BICYCLE_BINARY_DIR}
                        COMMENT "Generating API Documentation with Doxygen"
                        VERBATIM)
    else()
      message(FATAL_ERROR "Doxygen executable not found.")
    endif()
endif()

if (BICYCLE_BUILD_PYTHON_EXTENSION_MODULES)
  find_package(PythonLibs 3.2 REQUIRED)  # I am using the 3.2 C API
  find_package(PythonInterp 3.2 REQUIRED)# I am using the 3.2 C API
  find_package(NumPy 1.8 REQUIRED)       # I am using the new NumPy API
  add_subdirectory(python_wrapper)
endif()
