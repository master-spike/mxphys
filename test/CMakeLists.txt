cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0069 NEW) 

project(MxPhysTest)

add_executable(collision_conservation "${CMAKE_CURRENT_SOURCE_DIR}/collision_conservation.cpp")
set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED 20)
set_property(GLOBAL PROPERTY CXX_STANDARD 20)
set_property(TARGET mxphysdemo PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)



include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if(result)
    set_property(TARGET collision_conservation PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

target_include_directories(collision_conservation PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_link_libraries(collision_conservation PUBLIC MxPhysLib)
