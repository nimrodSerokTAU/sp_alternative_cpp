cmake_minimum_required(VERSION 3.14)
project(sp_alternative CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/main.cpp
    src/utils.cpp
    src/node.cpp
    src/unrooted_tree.cpp
    src/rooted_tree.cpp
    src/neighbor_joining.cpp
    src/global_alignment.cpp
    src/sp_score.cpp
    src/distance_calc.cpp
    src/sop_stats.cpp
    src/w_sop_stats.cpp
    src/gaps_stats.cpp
    src/entropy_stats.cpp
    src/kmer_stats.cpp
    src/tree_stats.cpp
    src/dist_labels_stats.cpp
    src/msa.cpp
    src/multi_msa_service.cpp
)

add_executable(sp_alternative ${SOURCES})

# Link filesystem library (needed on some Linux systems with older GCC)
if(NOT APPLE)
    target_link_libraries(sp_alternative PRIVATE stdc++fs)
endif()
