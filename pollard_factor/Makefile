CXX=mpic++
CXXFLAGS=-O3 -I../common
all: n_gen pollard_factor_test
n_gen: n_gen.o
	mpic++ -O3 -o n_gen n_gen.o -lntl -lgmp
common:
	make -C ../common all
pollard_factor_test: common pollard_factor_test.o
	mpic++ -o pollard_factor_test pollard_factor_test.o ../common/pollard_factor.o -lntl -lgmp
clean:
	rm -rf *.o pollard_factor_test n_gen *~
