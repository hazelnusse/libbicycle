project(BICYCLE CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${BICYCLE_SOURCE_DIR}/cmake/Modules)

option(BUILD_SHARED_LIBS "Build shared libraries." ON)
option(BICYCLE_BUILD_TESTS "Build tests." ON)
option(BICYCLE_BUILD_EXAMPLES "Build examples." ON)
option(BICYCLE_BUILD_DOXYGEN_DOCS "Build doxygen documentation." OFF)
option(BICYCLE_DOXYGEN_MATHJAX_USECDN "Use cdn servers to host MathJax (overrides BICYCLE_DOXYGEN_MATHJAX_LOCAL" ON)
set(BICYCLE_DOXYGEN_MATHJAX_LOCAL "" CACHE PATH "Location of MathJax javascript library")
option(BICYCLE_BUILD_PYTHON_EXTENSION_MODULES "Build extension modules for Python" OFF)

# Set g++ compiler flags if GNU C++ compiler is being used
if(CMAKE_COMPILER_IS_GNUCXX)
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if(GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
        set(ENABLE_CXX11 "-std=c++11")
    elseif(GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
        set(ENABLE_CXX11 "-std=c++0x")
    else()
        message(FATAL_ERROR "C++ compiler does not support C++11")
    endif()
    set(CMAKE_CXX_FLAGS "-Wall -Wno-unused-local-typedefs ${ENABLE_CXX11}")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
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

# If we find git installed, set BICYCLE_VERSION_COMMIT to most recent SHA1
find_package(Git)
if(GIT_FOUND)
  if(EXISTS ${BICYCLE_SOURCE_DIR}/.git)
    if(IS_DIRECTORY ${BICYCLE_SOURCE_DIR}/.git)
      execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --pretty=oneline
                    WORKING_DIRECTORY ${BICYCLE_SOURCE_DIR}
                    OUTPUT_VARIABLE BICYCLE_VERSION_COMMIT)

      # Keep everything up to the first space, add quotes
      string(REGEX REPLACE "([^ ]+).*"
             "\"\\1\"" BICYCLE_VERSION_COMMIT
             ${BICYCLE_VERSION_COMMIT})
    endif()
  endif()
else()
  set(BICYCLE_VERSION_COMMIT "\"\"")
endif()

if (BICYCLE_BUILD_DOXYGEN_DOCS)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
      add_subdirectory(doc)
    else()
      message(FATAL_ERROR "Doxygen executable not found.")
    endif()
endif()

if (BICYCLE_BUILD_PYTHON_EXTENSION_MODULES)
  find_package(PythonLibs 3.2 REQUIRED)  # I am using the 3.2 C API
  find_package(PythonInterp 3.2 REQUIRED)# I am using the 3.2 C API
  find_package(NumPy 1.8 REQUIRED)       # I am using the new NumPy API
  # Set a few variables used by the FindBoost.cmake script
  set(Boost_USE_MULTITHREADED     OFF)
  set(Boost_USE_STATIC_RUNTIME    OFF)
  find_package(Boost 1.52.0 COMPONENTS python)
  add_definitions(-DWRAP_PYTHON)
  add_subdirectory(python_wrapper)
endif()
