cmake_minimum_required(VERSION 2.8.3)
set(CMAKE_CXX_STANDARD 14)
project(libtetrahedra)

# set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-deprecated-declarations")

set(BUILD_SHARED_LIBS ON)

file(GLOB INCLUDE_FILES "include/*.h")

find_package(Eigen3 REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(
    ${EIGEN3_INCLUDE_DIRS} 
    ${PYTHON_INCLUDE_DIRS}
    matplotlib-cpp
    include)

add_executable(simple src/simple_construction.cpp)
target_link_libraries(simple ${PYTHON_LIBRARIES})