cmake_minimum_required(VERSION 3.19)
project(SPST_new)

set(CMAKE_CXX_STANDARD 14)

add_executable(main
        index/avl_tree.cpp index/spst.cpp cracking/cracking_util.cpp
        cracking/standard_cracking.cpp
        cracking/updates.cpp main.cpp index/skiplist_with_rank.h include/basestructs.h)