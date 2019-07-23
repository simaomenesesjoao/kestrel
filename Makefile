CC=g++ -std=c++11
CFLAGS=-I/usr/include/eigen3 -I${HOME}/include/eigen3
DEPS=src/*.hpp
OBJ=src/kpm_vector.o src/main.o src/aux.o src/hamiltonian.o src/cheb.o

%.o: %.cpp $(DEPS)
	$(CC) -Wall -DVERBOSE=$(VERBOSE) -DSTRIDE=$(STRIDE) -Dnis_complex=1 -O3 -c -o $@ $< $(CFLAGS)

kestrel: $(OBJ) 
	$(CC) -Wall -DVERBOSE=$(VERBOSE) -DSTRIDE=$(STRIDE) -Dnis_complex=1 -o $@ $^ $(CFLAGS)
