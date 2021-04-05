CC=g++ -fopenmp -pthread -std=c++11 -Wall
CFLAGS=-I/usr/include/eigen3 -I${HOME}/include/hdf5 -I${HOME}/include/eigen3 -I/usr/include -I${HOME}/include -I${HOME}/hdf_head
DEPS=src/*.hpp
OBJ=src/kpm_vector.o src/main.o src/aux.o src/hamiltonian.o src/cheb.o src/ComplexTraits.o src/myHDF5.o
LIB=-L/opt/local/lib -L/${HOME}/hdf_lib_bin  -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5

%.o: %.cpp $(DEPS)
	$(CC) -DMODEL=$(MODEL) -DVERBOSE=$(VERBOSE) -DSTRIDE=$(STRIDE) -Dnis_complex=1 -O3 -c -o $@ $< $(CFLAGS)

kestrel: $(OBJ) 
	$(CC) -DMODEL=$(MODEL) -DVERBOSE=$(VERBOSE) -DSTRIDE=$(STRIDE) -Dnis_complex=1 $(LIB) -o $@ $^ $(CFLAGS)
