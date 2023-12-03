CXX=g++
CXXFLAGS=-I
DEPS = generals.h lsh_func.h
OBJ1 = lsh.o generals.o lsh_func.o
OBJ2 = cube.o generals.o lsh_func.o cube_func.o
OBJ3 = cluster.o cluster_func.o generals.o lsh_func.o cube_func.o
OBJ4 = graph_search.o graph_search_func.o generals.o lsh_func.o 
%.o: %.c ${DEPS}
	$(CXX) -c -o $@ $< $(CXXFLAGS)

lsh: $(OBJ1)
	$(CXX) -o $@ $(OBJ1)

cube: $(OBJ2)
	$(CXX) -o $@ $(OBJ2)

cluster: $(OBJ3)
	$(CXX) -o $@ $(OBJ3)

graph_search: $(OBJ4)
	$(CXX) -pthread -o $@ $(OBJ4)

lsh.o : lsh.cpp 
	${CXX} -c lsh.cpp

cube.o : cube.cpp
	${CXX} -c cube.cpp

cluster.o : cluster.cpp
	${CXX} -c cluster.cpp

graph_search.o : graph_search.cpp
	${CXX} -c graph_search.cpp

graph_search_func.o : graph_search_func.cpp graph_search_func.h
	${CXX} -c graph_search_func.cpp graph_search_func.h

generals.o : generals.cpp generals.h
	${CXX} -c generals.cpp generals.h

lsh_func.o : lsh_func.cpp lsh_func.h
	${CXX} -c lsh_func.cpp lsh_func.h

cube_func.o : cube_func.cpp cube_func.h
	${CXX} -c cube_func.cpp cube_func.h

cluster_func.o : cluster_func.cpp cluster_func.h
	${CXX} -c cluster_func.cpp cluster_func.h

.PHONY: clean //necessary in case file with name clean exists
clean:
	rm -f *.o lsh cube
