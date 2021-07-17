COMPILER = g++
# *** object files ***
OBJECTS = utilityfunctions.o geometricfunctions.o enclosingtriangle.o delaunaytriangulation.o points.o main.o

FLAGS = -O2
# for debugging: -g
DEBUG = 

%.o : %.cpp
	$(COMPILER) -c $< $(FLAGS) $(DEBUG)

# *** TIN: Delaunay triangulation ***
TIN : $(OBJECTS)
	$(COMPILER) -o $@ $(OBJECTS) $(FLAGS)

clean:
	rm -f $(OBJECTS) TIN
