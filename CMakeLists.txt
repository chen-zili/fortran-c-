# 设置项目所需cmake的最低版本
cmake_minimum_required(VERSION 3.0.0)

# 项目名称，也是最后生成的可执行文件名
project(fortran-cpp)

# 启用对Fortran语言的支持
enable_language(Fortran C CXX)

# 设置全局路径
set (PROJECTPATH ${PROJECT_SOURCE_DIR})

# add_compile_options(-fPIC)

set(CMAKE_BUILD_TYPE Debug CACHE STRING "set build type to debug")

include_directories(${PROJECT_SOURCE_DIR}/transLine/include)
include_directories(${PROJECT_SOURCE_DIR}/petsc/include/)
include_directories(${PROJECT_SOURCE_DIR}/petsc/include-double/)
include_directories(${PROJECT_SOURCE_DIR}/fmt/include/)
include_directories(/opt/intel/oneapi/mpi/latest/include/)

file(GLOB_RECURSE SRCS ${PROJECT_SOURCE_DIR}/*.f90 tem1.c)

# 添加源文件
add_executable(${PROJECT_NAME} ${SRCS})

add_subdirectory(transLine)

# 链接
target_link_libraries(${PROJECT_NAME} transLine)
