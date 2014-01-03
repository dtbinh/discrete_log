CXX=mpic++
CXXFLAGS=-O3
LIBS=-lntl -lgmp
all: common

clean:
	rm -rf *~ *.o

common:	gauss.o mpi_ntl_routines.o sieve.o pollard_factor.o

mpi_ntl_routines_test:
	$(CXX) $(CXXFLAGS) -DMODULE_TEST -o mpi_ntl_routines mpi_ntl_routines.cpp $(LIBS)
