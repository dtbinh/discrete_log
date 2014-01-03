CXX=mpic++
CXXFLAGS=-O3 -I../common
all: index_calc task_gen
index_calc: common index_calc.o
	mpic++ -O3 -o index_calc index_calc.o ../common/gauss.o ../common/mpi_ntl_routines.o ../common/sieve.o -lntl -lgmp
common:
	make -C ../common all
task_gen:task_gen.o
	g++ -o task_gen task_gen.o -lntl -lgmp
clean:
	rm -rf *.o index_calc task_gen *~
