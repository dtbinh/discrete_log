#ifndef NONSMOOTH_EXT_H
#define NONSMOOTH_EXT_H 1
#include <NTL/ZZ.h>
#include "mpi.h"
#include "mpi_ntl_routines.h"
NTL_CLIENT

bool make_smooth(ZZ x, ZZ p, ZZ smooth_checker, ZZ * basis, long basis_size, //in
		 int crank, int csize, int block_size,  //in
		ZZ& rez, ZZ& rez_power, ZZ * result_powers); //out


bool solve_matrix(ZZ ** base, int lines, int columns, ZZ r, //in
		  ZZ& alpha, ZZ& beta); //out

#endif
