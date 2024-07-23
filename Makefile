SRC = src
BIN = bin

# List sources manually to exclude sim.cpp and simmulti.cpp
SOURCE = decay.cpp frag.cpp Gobbiarray.cpp loss.cpp polyScat.cpp mScat.cpp random.cpp frame.cpp tele.cpp correlations.cpp rootoutput.cpp
OBJECT = $(patsubst %, $(BIN)/%, $(notdir $(SOURCE:.cpp=.o)))
HEADER = rootoutput.h correlations.h
DICTSRC = $(patsubst %, $(SRC)/%, $(HEADER))
LIB = /home/Li6Webb/simlib/simlib.a

CC = c++
LINKOPTION = $(shell root-config --libs --cflags)
INCLUDE = -I$(shell root-config --incdir) -I/home/Li6Webb/simlib
CFLAGS = -c -w -std=c++14 $(INCLUDE)

sim : $(BIN)/sim_Li6_alphapn.o $(OBJECT) libLi6sim.cxx
	@echo "Linking..."
	$(CC) -o$@ $(filter-out %.cxx, $^) $(LIB) $(LINKOPTION) $(filter %.cxx, $^)
	@echo "Finished"

sim_3+ : $(BIN)/sim_Li6_dalpha.o $(OBJECT) libLi6sim.cxx
	@echo "Linking..."
	$(CC) -o$@ $(filter-out %.cxx, $^) $(LIB) $(LINKOPTION) $(filter %.cxx, $^)
	@echo "Finished"

$(BIN)/%.o : $(SRC)/%.cpp
	@echo "Compiling..."
	$(CC) $(CFLAGS) $< -o$@

libLi6sim.cxx : $(DICTSRC) $(SRC)/LinkDef.h
	@echo "Creating dictionary..."
	rootcling $@ $^

clean:
	rm -f $(BIN)/*.o
	rm -f *.cxx
	rm -f *.pcm
	rm -f sim
	rm -f sim_3+
	
help:
	@echo "commands ->"
	@echo "make"
	@echo "  - equivalent to \"make sim\""
	@echo "make sim_3+"
	@echo "make clean"
