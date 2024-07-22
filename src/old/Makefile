
SRC = src
BIN = bin

SOURCE = $(widcard $(SRC)/*.cpp)
OBJECT = $(patsubst %, $(BIN)/%, $(notdir $(SOURCE:.cpp=.o)))
#OBJECTS = decay.o frag.o fragneut.o Gobiarray.o Neutarray.o loss.o mScat.o random.o frame.o fragment.o tele.o correlations.o

#ALLOBJECTS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(wildcard *.cpp))
#DEPS := $(patsubst %.cpp,%.h,$(wildcard *.cpp))
#CFLAGS= -c -std=c++11 -I$(shell root-config --incdir)
#COMPILER= c++
#LINKOPTION = $(shell root-config --libs) 

sim: sim.o $(OBJDIR)/$(OBJECTS)
	$(COMPILER) -o sim sim.o $(OBJDIR)/$(OBJECTS) $(LINKOPTION)
	
multi: simmulti.o $(OBJDIR)/$(OBJECTS)
	$(COMPILER) -o simmulti simmulti.o $(OBJECTS) $(LINKOPTION)
 
$(OBJDIR)/$(ALLOBJECTS): %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@
	

clean:
	rm -f *.o *.*~
	
help:
	@echo "src: $(SOURCE)"
	@echo "obj: $(OBJECT)"
