# 设置项目所需cmake的最低版本
cmake_minimum_required(VERSION 3.0.0)

# 项目名称，也是最后生成的可执行文件名
project(fortran-cpp)

# 启用对Fortran语言的支持
enable_language(Fortran C CXX)
set(CMAKE_CXX_STANDARD 11)

# add_compile_options(-fPIC)

set(CMAKE_BUILD_TYPE Debug CACHE STRING "set build type to debug")

# 设置源文件和头文件路径
include_directories(${PROJECT_SOURCE_DIR}/)
file(GLOB_RECURSE SRCS ${PROJECT_SOURCE_DIR}/*.f90 ${PROJECT_SOURCE_DIR}/*.cpp)

# 添加源文件
add_executable(${PROJECT_NAME} ${SRCS})

