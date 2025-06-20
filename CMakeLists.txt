cmake_minimum_required(VERSION 3.16)
project(li6sim_wrapper)

# Option to build full project or just the library
option(LI6SIM_LIBRARY_ONLY "Only build the LI6SIM static library" OFF)

if(LI6SIM_LIBRARY_ONLY)
    include(${CMAKE_CURRENT_LIST_DIR}/CMakeListsLib.txt)
else()
    include(${CMAKE_CURRENT_LIST_DIR}/CMakeListsFull.txt)
endif()
