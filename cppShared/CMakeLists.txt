cmake_minimum_required(VERSION 3.0.0)

project(cppShared)

# 支持c++ 11
enable_language(CXX C)
set(CMAKE_CXX_STANDARD 11)

# 设置头文件路径
include_directories(${CMAKE_CURRENT_LIST_DIR}/include/)
file(GLOB_RECURSE LIB_SRC ${CMAKE_CURRENT_LIST_DIR}/source/*.cpp)

# 生成
add_library(cppShared SHARED ${LIB_SRC})

set_target_properties(cppShared PROPERTIES OUTPUT_NAME "cppShared")

set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_LIST_DIR}/lib/)

