cmake_minimum_required(VERSION 2.8)
PROJECT(maps)
set(exe_name maps)

FIND_PACKAGE(OpenCV REQUIRED)

include_directories(.)

file(GLOB cpp_files *.cpp)
file(GLOB h_files *.h)
add_executable(${exe_name} ${cpp_files} ${h_files})
target_link_libraries(${exe_name} ${OpenCV_LIBRARIES})