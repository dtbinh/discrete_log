#ifndef MPI_EC_ROUTINES_H
#define MPI_EC_ROUTINES_H 1

#include "mpi_ntl_routines.h"
#include "ec_point.h"
#include "mpi.h"
#include <NTL/ZZ.h>

NTL_CLIENT

//Послать EC точку
int MPI_EC_Send(ec_point  x, int process_id);
//Принять EC точку
int MPI_EC_Recv(ec_point& x, int process_id);
//Широковещательно разослать EC точку
int MPI_EC_Bcast(ec_point& x, int process_id);

void EC_Pack(ec_point x, unsigned char * buf, int& buf_ub);

bool EC_Unpack(ec_point& x, unsigned char * buf, int& buf_ub, int buf_size);

#endif
