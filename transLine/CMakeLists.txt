# 设置项目所需cmake的最低版本
cmake_minimum_required(VERSION 3.0.0)

# 项目名称，也是最后生成的可执行文件名
project(transLine)

# 语言支持
enable_language(CXX C)

# 支持c++ 11
set(CMAKE_CXX_STANDARD 11)

# 设置源文件和头文件路径
file(GLOB_RECURSE SRCS ${CMAKE_CURRENT_LIST_DIR}/source/*.cpp)

# 设置头文件路径
include_directories(${CMAKE_CURRENT_LIST_DIR}/include/)
include_directories(${PROJECTPATH}/petsc/include/)
include_directories(${PROJECTPATH}/petsc/include-double/)
include_directories(${PROJECTPATH}/fmt/include/)
include_directories(/opt/intel/oneapi/mpi/latest/include/)

# 生成
add_library(transLine SHARED ${SRCS})

# 链接库
target_link_libraries(transLine ${PROJECTPATH}/petsc/lib/libpetsc.so)
target_link_libraries(transLine ${PROJECTPATH}/fmt/lib/libfmt.so.8)
target_link_libraries(transLine /opt/intel/oneapi/mpi/2021.3.0/lib/release/libmpi.so.12)
