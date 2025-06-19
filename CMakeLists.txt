cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(li6sim)

# Setup directories.
set(SRC ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_INSTALL_LIBDIR ${PROJECT_BINARY_DIR}/lib)
set(SIMLIB /home/Li6Webb/simlib)

# Set project sources
set(SOURCES decay.cpp frag.cpp Gobbiarray.cpp loss.cpp polyScat.cpp mScat.cpp random.cpp frame.cpp tele.cpp correlations.cpp rootoutput.cpp)
set(LIBHEADERS rootoutput.h correlations.h)

list(TRANSFORM SOURCES PREPEND ${SRC}/)
list(TRANSFORM LIBHEADERS PREPEND ${SRC}/)

# Locate the ROOT package and define a number of useful targets and variables
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Define cross-section file path
set(XSECPATH ${CMAKE_CURRENT_SOURCE_DIR}/input)
add_definitions(-DXSECPATH=\"${XSECPATH}\")

# Create executable
add_executable(sim ${SRC}/sim_Li6_alphapn.cpp ${SOURCES})
target_link_libraries(sim ${SIMLIB}/simlib.a ROOT::RIO ROOT::Tree ROOT::Hist)
target_include_directories(sim PRIVATE ${SIMLIB})
set_target_properties(sim PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Generate dictionaries and add them as a dependency of the executable (via the MODULE parameter)
ROOT_GENERATE_DICTIONARY(libLi6sim ${LIBHEADERS} MODULE sim LINKDEF ${SRC}/LinkDef.h)

# Create executable
add_executable(sim_3+ ${SRC}/sim_Li6_dalpha.cpp ${SOURCES} ${PROJECT_SOURCE_DIR}/libLi6sim.cxx)
target_link_libraries(sim_3+ ${SIMLIB}/simlib.a ROOT::RIO ROOT::Tree ROOT::Hist)
target_include_directories(sim_3+ PRIVATE ${SIMLIB})
set_target_properties(sim_3+ PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
