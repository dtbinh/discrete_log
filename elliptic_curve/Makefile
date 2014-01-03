CXX=mpic++
CXXFLAGS=-O3 -I../common
LIBS=-lntl -lgmp
all: ec_pollard ec_task_gen meet_in_the_middle meet_on_random_tree
common:
	make -C ../common
ec_pollard: ec_point.o ec_routine.o ec_pollard.o
	$(CXX) $(CXXFLAGS) -o ec_pollard ec_point.o ec_routine.o ec_pollard.o $(LIBS)
ec_task_gen: common ec_point.o ec_routine.o ec_task_gen.o
	$(CXX) $(CXXFLAGS) -o ec_task_gen ec_point.o ec_routine.o ec_task_gen.o ../common/pollard_factor.o $(LIBS)
meet_in_the_middle: ec_point.o ec_routine.o meet_in_the_middle.o mpi_ec_routines.o common
	$(CXX) $(CXXFLAGS) -o meet_in_the_middle meet_in_the_middle.o ec_point.o ec_routine.o mpi_ec_routines.o ../common/mpi_ntl_routines.o $(LIBS)
meet_on_random_tree: ec_point.o ec_routine.o meet_on_random_tree.o mpi_ec_routines.o rb_tree.o common
	$(CXX) $(CXXFLAGS) -o meet_on_random_tree meet_on_random_tree.o rb_tree.o ec_point.o ec_routine.o mpi_ec_routines.o ../common/mpi_ntl_routines.o $(LIBS)
<<<<<<< .mine
#mitm_tune_param: mitm_tune_param.o
#	$(CXX) $(CXXFLAGS) -o mitm_tune_param mitm_tune_param.o $(LIBS)
=======
mort2: ec_point.o ec_routine.o mort2.o mpi_ec_routines.o rb_tree.o common
	$(CXX) $(CXXFLAGS) -o mort2 mort2.o rb_tree.o ec_point.o ec_routine.o mpi_ec_routines.o ../common/mpi_ntl_routines.o $(LIBS)
	
mitm_tune_param: mitm_tune_param.o
	$(CXX) $(CXXFLAGS) -o mitm_tune_param mitm_tune_param.o $(LIBS)
>>>>>>> .r1466
mpi_ec_test: common ec_point.o ec_routine.o mpi_ec_routines.o mpi_ec_test.o
	$(CXX) $(CXXFLAGS) -o mpi_ec_test ec_point.o ec_routine.o mpi_ec_routines.o mpi_ec_test.o ../common/mpi_ntl_routines.o $(LIBS)
cup_table_test: cup_table.o cup_table_test.o
	$(CXX) $(CXXFLAGS) -o cup_table_test cup_table.o cup_table_test.o $(LIBS)
clean:
	rm -rf *.o ec_pollard ec_task_gen mpi_ec_test meet_in_the_middle meet_on_random_tree *~
