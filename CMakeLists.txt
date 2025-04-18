cmake_minimum_required(VERSION 3.24)
project(Jarrow)

set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")

include_directories("/opt/homebrew/opt/llvm/include")
link_directories("/opt/homebrew/opt/llvm/lib")

# Create executable target first
add_executable(Jarrow main.cpp)

# Then link OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Jarrow PUBLIC OpenMP::OpenMP_CXX)
endif()