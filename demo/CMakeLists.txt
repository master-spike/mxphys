cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0069 NEW) 

project(MxPhysDemo)

find_package(SDL2 REQUIRED)

add_executable(mxphysdemo "${CMAKE_CURRENT_SOURCE_DIR}/demo.cpp")
add_executable(mxphysdemo_rope "${CMAKE_CURRENT_SOURCE_DIR}/demo_ropes.cpp")
set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED 20)
set_property(GLOBAL PROPERTY CXX_STANDARD 20)

include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if(result)
    set_property(TARGET mxphysdemo PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET mxphysdemo_rope PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

target_include_directories(mxphysdemo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_link_libraries(mxphysdemo PUBLIC SDL2::SDL2 SDL2::SDL2main MxPhysLib)


target_include_directories(mxphysdemo_rope PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_link_libraries(mxphysdemo_rope PUBLIC SDL2::SDL2 SDL2::SDL2main MxPhysLib)
