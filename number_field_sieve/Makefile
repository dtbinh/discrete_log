CXX=mpic++
CXXFLAGS=-O3 -I../common
LIBS=-lntl -lgmp
all: number_field_sieve smooth_task_gen
number_field_sieve: common number_field_sieve.o algebraic_number.o nonsmooth_ext.o
	$(CXX) -O3 -o number_field_sieve number_field_sieve.o nonsmooth_ext.o ../common/mpi_ntl_routines.o algebraic_number.o ../common/gauss.o ../common/sieve.o $(LIBS)
common:
	make -C ../common all
smooth_task_gen: smooth_task_gen.o
	$(CXX) -O3 -o smooth_task_gen smooth_task_gen.o $(LIBS)
nonsmooth_ext:common
	$(CXX) -I../common -DMODULE_TEST -O3 -o nonsmooth_ext nonsmooth_ext.cpp ../common/sieve.o $(LIBS)
clean:
	rm -rf *.o number_field_sieve smooth_task_gen *~
