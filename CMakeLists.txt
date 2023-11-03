cmake_minimum_required(VERSION 3.8)

project(MxPhysProject)

if (MSVC)
    SET(CMAKE_CXX_FLAGS "/W3 /WX /O2 /EHsc")
else()
    SET(CMAKE_CXX_FLAGS "-Wall -Wextra -Wconversion -Wpedantic -Werror -O3")
endif()

add_subdirectory(src)
add_subdirectory(demo)
add_subdirectory(test)
