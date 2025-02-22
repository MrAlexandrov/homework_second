cmake_minimum_required(VERSION 3.14)
project(homework_second CXX)

set(CMAKE_CXX_STANDARD 20)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

file(GLOB SRC_FILES "src/*.cpp")

add_executable(${PROJECT_NAME} main.cpp ${SRC_FILES})

# enable_testing()
# add_subdirectory(${CMAKE_SOURCE_DIR}/tests)

# set(CTEST_BINARY_DIRECTORY ${CMAKE_BINARY_DIR})
