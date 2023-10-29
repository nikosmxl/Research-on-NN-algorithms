

CXX=g++
CXXFLAGS=-I
DEPS = generals.h lsh_func.h
OBJ1 = lsh.o generals.o lsh_func.o
OBJ2 = cube.o generals.o lsh_func.o
%.ο: %.c ${DEPS}
	$(CC) –c –o $@ $< $(CFLAGS)
lsh: $(OBJ1)
	$(CC) –o $@ $^ $(CFLAGS)

cube: $(OBJ2)
	$(CC) –o $@ $^ $(CFLAGS)

.PHONY: clean //necessary in case file with name clean exists
clean:
	rm –f *.o lsh cube