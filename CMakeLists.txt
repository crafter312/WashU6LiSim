cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(li6sim)

# Option to build full project or just the library
option(LI6SIM_LIBRARY_ONLY "Only build the LI6SIM static library" OFF)

# Enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Setup directories.
set(SRC ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_INSTALL_LIBDIR ${PROJECT_BINARY_DIR}/lib)

set(USER_HOME_DIRECTORY $ENV{HOME})
message(STATUS "User Home Directory: ${USER_HOME_DIRECTORY}")
set(SIMLIB ${USER_HOME_DIRECTORY}/simlib)

# Set project sources
set(SOURCES decay.cpp frag.cpp Gobbiarray.cpp loss.cpp polyScat.cpp mScat.cpp random.cpp frame.cpp tele.cpp correlations.cpp rootoutput.cpp)
set(LIBHEADERS rootoutput.h correlations.h)

list(TRANSFORM SOURCES PREPEND ${SRC}/)
list(TRANSFORM LIBHEADERS PREPEND ${SRC}/)

# Locate the ROOT package and define a number of useful targets and variables
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist Graf)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

# Define cross-section file path
set(XSECPATH ${CMAKE_CURRENT_SOURCE_DIR}/input)
add_definitions(-DXSECPATH=\"${XSECPATH}\")

# Define loss file path
set(LOSSPATH ${CMAKE_CURRENT_SOURCE_DIR}/loss)
add_definitions(-DLOSSPATH=\"${LOSSPATH}\")

# Generate the dictionary using ROOT_GENERATE_DICTIONARY
ROOT_GENERATE_DICTIONARY(G__LI6SIM ${LIBHEADERS} MODULE LI6SIM LINKDEF ${SRC}/LinkDef.h)

# Create a custom target for the dictionary generation
add_custom_target(generate_dict ALL DEPENDS ${CMAKE_BINARY_DIR}/G__LI6SIM.cxx)

# Create an OBJECT library to handle the generated sources (with sim class)
add_library(Li6simObjects OBJECT ${SRC}/Li6sim_alphapn.cpp ${SOURCES} ${CMAKE_BINARY_DIR}/G__LI6SIM.cxx)
set_target_properties(Li6simObjects PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(Li6simObjects PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(Li6simObjects ${SIMLIB}/simlib.a ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf)
target_include_directories(Li6simObjects PRIVATE ${SIMLIB})
add_dependencies(Li6simObjects generate_dict)

# Create library
add_library(LI6SIM SHARED $<TARGET_OBJECTS:Li6simObjects>)
set_target_properties(LI6SIM PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(LI6SIM PROPERTIES OUTPUT_NAME "LI6SIM")
target_link_libraries(LI6SIM ${SIMLIB}/simlib.a ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf)
target_include_directories(LI6SIM PUBLIC ${SIMLIB})
target_include_directories(LI6SIM PUBLIC ${SRC} ${ROOT_INCLUDE_DIRS})
add_dependencies(LI6SIM generate_dict)

if(LI6SIM_LIBRARY_ONLY)
else()
	# Define dictionary file path
	set(PCMFILE ${CMAKE_BINARY_DIR}/libLI6SIM_rdict.pcm)
	add_definitions(-DPCMFILE=\"${PCMFILE}\")

	# Define shared library file path (this also forces the correct file name)
	set(SOFILE ${CMAKE_BINARY_DIR}/libLI6SIM.so)
	add_definitions(-DSOFILE=\"${SOFILE}\")

	# Create 6Li(2+2)->alpha+p+n executable
	add_executable(sim ${SRC}/sim_Li6_alphapn.cpp $<TARGET_OBJECTS:Li6simObjects>)
	target_link_libraries(sim ${SIMLIB}/simlib.a ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf)
	target_include_directories(sim PRIVATE ${SIMLIB})
	set_target_properties(sim PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
	add_dependencies(sim generate_dict)

	# Create 6Li(3+)->alpha+d executable
	add_executable(sim_3+ ${SRC}/sim_Li6_dalpha.cpp $<TARGET_OBJECTS:Li6simObjects>)
	target_link_libraries(sim_3+ ${SIMLIB}/simlib.a ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf)
	target_include_directories(sim_3+ PRIVATE ${SIMLIB})
	set_target_properties(sim_3+ PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
	add_dependencies(sim_3+ generate_dict)
endif()
