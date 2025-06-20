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

# Generate ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__Li6sim ${LIBHEADERS} MODULE Li6sim LINKDEF ${SRC}/LinkDef.h)

# Create a custom target to generate the dictionary
add_custom_target(generate_dict ALL
    COMMAND ${CMAKE_COMMAND} -E echo "Generating libLi6sim dictionary"
    DEPENDS ${CMAKE_BINARY_DIR}/G__Li6sim.cxx
)

# Create 6Li(2+2)->alpha+p+n executable
add_executable(sim ${SRC}/sim_Li6_alphapn.cpp ${SRC}/Li6sim_alphapn.cpp ${SOURCES} ${PROJECT_SOURCE_DIR}/G__Li6sim.cxx)
target_link_libraries(sim ${SIMLIB}/simlib.a ROOT::RIO ROOT::Tree ROOT::Hist)
target_include_directories(sim PRIVATE ${SIMLIB})
set_target_properties(sim PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_dependencies(sim generate_dict)

# Create 6Li(3+)->alpha+d executable
add_executable(sim_3+ ${SRC}/sim_Li6_dalpha.cpp ${SOURCES} ${PROJECT_SOURCE_DIR}/G__Li6sim.cxx)
target_link_libraries(sim_3+ ${SIMLIB}/simlib.a ROOT::RIO ROOT::Tree ROOT::Hist)
target_include_directories(sim_3+ PRIVATE ${SIMLIB})
set_target_properties(sim_3+ PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_dependencies(sim_3+ generate_dict)

# Create library
add_library(LI6SIM STATIC ${SRC}/Li6sim_alphapn.cpp ${SOURCES} ${PROJECT_SOURCE_DIR}/G__Li6sim.cxx)
set_target_properties(LI6SIM PROPERTIES EXCLUDE_FROM_ALL TRUE)
target_link_libraries(LI6SIM ${SIMLIB}/simlib.a ROOT::RIO ROOT::Tree ROOT::Hist)
target_include_directories(LI6SIM PRIVATE ${SIMLIB})
target_include_directories(LI6SIM PUBLIC ${SRC} ${ROOT_INCLUDE_DIRS})
add_dependencies(LI6SIM generate_dict)
