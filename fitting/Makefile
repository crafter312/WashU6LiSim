OBJECTS = fit.o
ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
CFLAGS= -c -ggdb -std=c++14 -I$(shell root-config --incdir)  -I/usr/local/include/gsl $(shell gsl-config --cflags) -I/home/Li7star/simlib 
COMPILER= c++
LINKOPTION = -lgsl -lgslcblas $(shell root-config --libs)  
LIB = /home/Li7star/simlib/simlib.a

chi: chi.o $(OBJECTS)
	$(COMPILER) -o chi chi.o $(OBJECTS) $(LIB) $(LINKOPTION)

chi_multi: chi_multi.o $(OBJECTS)
	$(COMPILER) -o chi_multi chi_multi.o $(OBJECTS) $(LIB) $(LINKOPTION)
 
$(ALLOBJECTS): %.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@


clean:
	rm -f *.o

