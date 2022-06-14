cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKS_CXX_STANDARD_REQUIRED ON)

project(main VERSION 1.0)
file(GLOB_RECURSE SRC_FILES ../../src/FE8nodeICM.cpp)
add_executable(main main.cpp ${SRC_FILES})
target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
