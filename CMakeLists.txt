cmake_minimum_required(VERSION 3.16)
project(IIT)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-g -O3 -march=native -Wall -fopenmp -mavx -pthread -lrocksdb)

add_subdirectory(graph)
add_subdirectory(query)
add_subdirectory(utility)
#add_subdirectory(utils)
add_subdirectory(roaring)