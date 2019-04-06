cmake_minimum_required(VERSION 3.5)
project(mpm)

set(CMAKE_CXX_STANDARD 11)

include_directories(external/Eigen/)
include_directories(include)

file(GLOB HEADER_FILES include/*.h)
file(GLOB SOURCE_FILES src/*.cpp)

add_executable(mpm ${HEADER_FILES} ${SOURCE_FILES})