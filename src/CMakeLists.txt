cmake_minimum_required(VERSION 3.8)

include_directories("$(CMAKE_CURRENT_SOURCE_DIR)/../include")

add_library(MxPhysLib
    forces.cpp
    polygon.cpp
    body.cpp
    eventmanager.cpp
)

target_compile_features(MxPhysLib PUBLIC cxx_std_20)

set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED 20)
set_property(GLOBAL PROPERTY CXX_STANDARD 20)

SET(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -O3")
