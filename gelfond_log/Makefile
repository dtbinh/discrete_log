CXX=mpic++
CXXFLAGS=-O3 -I../common
all: gelfond_make_task gelfond_solve_task
gelfond_make_task: gelfond_make_task.o common
	mpic++ -O3 -o gelfond_make_task gelfond_make_task.o ../common/pollard_factor.o -lntl -lgmp
gelfond_solve_task: gelfond_solve_task.o
	mpic++ -O3 -o gelfond_solve_task gelfond_solve_task.o -lntl -lgmp
common:
	make -C ../common all
clean:
	rm -rf *.o gelfond_make_task gelfond_solve_task *~
