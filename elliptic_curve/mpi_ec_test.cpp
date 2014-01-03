#include "ec_point.h"
#include "ec_routine.h"
#include "pollard_factor.h"
#include <time.h>

#include <mpi.h>
#include "mpi_ec_routines.h"
#include <NTL/ZZ.h>
#include <iostream>
using namespace std;

int main(int argc, char ** argv) {
MPI_Init(&argc,&argv);
int crank,csize;
MPI_Comm_rank(MPI_COMM_WORLD,&crank);
MPI_Comm_size(MPI_COMM_WORLD,&csize);
//Ввести параметры кривой
ec_point::Init(to_ZZ(1),to_ZZ(1),to_ZZ(657469));
ec_point P(to_ZZ(72342),to_ZZ(448573)),Q(to_ZZ(608539),to_ZZ(217169)),inf(true);
ec_point p1,q1,inf1;

if (crank==1)
	cout << "Testing P2P connections - normal: ";
if (crank==0) {
	P.auxP=1;
	P.auxQ=2;
	MPI_EC_Send(P,1);
	}
else if (crank==1) {
	MPI_EC_Recv(p1,0);
	if (p1==P && p1.auxP==1 && p1.auxQ==2)
		cout << "pass" <<endl;
	else 
		cout << "FAIL" <<endl;
	}

MPI_Barrier(MPI_COMM_WORLD);

if (crank==0)
	cout << "Testing P2P connections - inf: ";
if (crank==1) {
	inf.auxP=10;
	inf.auxQ=20;
	MPI_EC_Send(inf,0);
	}
else if (crank==0) {
	MPI_EC_Recv(inf1,1);
	if (inf1==to_ec_point(true) && inf1.auxP==10 && inf1.auxQ==20)
		cout << "pass" <<endl;
	else 
		cout << "FAIL" <<endl;
	}

MPI_Barrier(MPI_COMM_WORLD);

if (crank==1)
	cout << "Testing Bcast connections - normal: ";
if (crank==0) {
	P=Q;
	P.auxP=11;
	P.auxQ=12;
	}

MPI_EC_Bcast(P,0);

if (crank==1) {
	if (P==Q && P.auxP==11 && P.auxQ==12)
		cout << "pass" <<endl;
	else 
		cout << "FAIL" <<endl;
	}

MPI_Barrier(MPI_COMM_WORLD);

if (crank==0)
	cout << "Testing Bcast connections - inf: ";

if (crank==1) {
	P=inf;
	P.auxP=110;
	P.auxQ=120;
	}
MPI_EC_Bcast(P,1);

if (crank==0) {
	if (P==to_ec_point(true) && P.auxP==110 && P.auxQ==120)
		cout << "pass" <<endl;
	else 
		cout << "FAIL" <<endl;
	}

MPI_Barrier(MPI_COMM_WORLD);


MPI_Finalize();
return 0;
}
