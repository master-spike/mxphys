cmake_minimum_required(VERSION 3.5)

project(MxPhysDemo)

add_executable(demo "${CMAKE_CURRENT_SOURCE_DIR}/demo/demo.cpp")
set_property(TARGET demo PROPERTY CXX_STANDARD_REQUIRED 20)
set_property(TARGET demo PROPERTY CXX_STANDARD 20)

target_include_directories(demo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(demo SDL2)


