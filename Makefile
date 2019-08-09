CC=g++ -pthread -std=c++11 -Wall
CFLAGS=-I/usr/include/eigen3 -I${HOME}/include/eigen3 -I/usr/include
DEPS=src/*.hpp
OBJ=src/kpm_vector.o src/main.o src/aux.o src/hamiltonian.o src/cheb.o src/ComplexTraits.o src/myHDF5.o src/tcp_client.o
LIB=-L/opt/local/lib  -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5

%.o: %.cpp $(DEPS)
	$(CC) -DVERBOSE=$(VERBOSE) -DSTRIDE=$(STRIDE) -Dnis_complex=1 -O3 -c -o $@ $< $(CFLAGS)

kestrel: $(OBJ) 
	$(CC) -DVERBOSE=$(VERBOSE) -DSTRIDE=$(STRIDE) -Dnis_complex=1 $(LIB) -o $@ $^ $(CFLAGS)
