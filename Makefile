CXX=g++
CXXFLAGS=-I
DEPS = generals.h lsh_func.h
OBJ1 = lsh.o generals.o lsh_func.o
OBJ2 = cube.o generals.o lsh_func.o
%.o: %.c ${DEPS}
	$(CXX) -c -o $@ $< $(CXXFLAGS)

lsh: $(OBJ1)
	$(CXX) -o $@ $(OBJ1)

cube: $(OBJ2)
	$(CXX) -o $@ $(OBJ2)

lsh.o : lsh.cpp 
	${CXX} -c lsh.cpp

cube.o : cube.cpp
	${CXX} -c cube.cpp

generals.o : generals.cpp generals.h
	${CXX} -c generals.cpp generals.h

lsh_func.o : lsh_func.cpp lsh_func.h
	${CXX} -c lsh_func.cpp lsh_func.h

.PHONY: clean //necessary in case file with name clean exists
clean:
	rm -f *.o lsh cube
