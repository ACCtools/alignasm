cmake_minimum_required(VERSION 3.27)
include(FetchContent)

project(alignasm)

set(CMAKE_CXX_STANDARD 20)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fprofile-arcs -ftest-coverage)
    add_link_options(--coverage)
endif()

find_package(unordered_dense REQUIRED)
find_package(TBB REQUIRED)

FetchContent_Declare(
        csv
        GIT_REPOSITORY https://github.com/vincentlaucsb/csv-parser.git
        GIT_SHALLOW TRUE
        GIT_TAG 2.2.1
)

FetchContent_Declare(
        argparse
        GIT_REPOSITORY https://github.com/p-ranav/argparse.git
        GIT_TAG v2.9
)

FetchContent_MakeAvailable(argparse)

FetchContent_MakeAvailable(csv)

add_executable(alignasm
        src/alignasm.cpp
        src/paf_data.hpp
        src/k_shortest_walks.hpp
        src/k_weighted_bfs.hpp
        src/priority_queue_vector.hpp
        src/graph_operations.hpp
        src/paf_data.cpp
        src/leftist_heap.hpp
)

target_include_directories(alignasm PRIVATE csv)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_libraries(alignasm gcov)
endif()

target_link_libraries(alignasm csv unordered_dense::unordered_dense TBB::tbb argparse)