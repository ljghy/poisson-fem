cmake_minimum_required(VERSION 3.8)

project(poisson-fem)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

add_subdirectory(src)
