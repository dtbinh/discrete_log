#include "mpi_ntl_routines.h"
#include "mpi_ec_routines.h"
#include "mpi.h"
#include <NTL/ZZ.h>

NTL_CLIENT

//Послать EC точку
int MPI_EC_Send(ec_point  x, int process_id) {
if (x.is_infinity()) {
	MPI_NTL_Send(to_ZZ(-1),process_id);
	}
else {
	MPI_NTL_Send(x.get_x(),process_id);
	MPI_NTL_Send(x.get_y(),process_id);
	}
MPI_NTL_Send(x.auxP,process_id);
MPI_NTL_Send(x.auxQ,process_id);
}

//Принять EC точку
int MPI_EC_Recv(ec_point& x, int process_id) {
ZZ new_x,new_y;
MPI_NTL_Recv(new_x,process_id);
if (new_x==to_ZZ(-1)) {
	x=to_ec_point(true);
	}
else {
	MPI_NTL_Recv(new_y,process_id);
	x=to_ec_point(new_x,new_y);
	}
MPI_NTL_Recv(x.auxP,process_id);
MPI_NTL_Recv(x.auxQ,process_id);
}

//Широковещательно разослать EC точку
int MPI_EC_Bcast(ec_point& x, int process_id) {
	bool is_inf = x.is_infinity();
	ZZ auxP=x.auxP,auxQ=x.auxQ;
	MPI_NTL_Bcast(auxP,process_id);
	MPI_NTL_Bcast(auxQ,process_id);

	MPI_Bcast(&is_inf,1,MPI_BYTE,process_id,MPI_COMM_WORLD);
	if (is_inf) {
		x = to_ec_point(true);
		}
	else {
		ZZ new_x=x.get_x(),new_y=x.get_y();
		MPI_NTL_Bcast(new_x,process_id);
		MPI_NTL_Bcast(new_y,process_id);
		x = to_ec_point(new_x,new_y);
		}
	x.auxP=auxP;
	x.auxQ=auxQ;
	}

void EC_Pack(ec_point x, unsigned char * buf, int& buf_ub) {
  if (x.is_infinity()) {
    buf[buf_ub++]=0;
    }
  else {
    buf[buf_ub++]=1;
    NTL_Pack(x.get_x(),buf,buf_ub);
    NTL_Pack(x.get_y(),buf,buf_ub);
    NTL_Pack(x.auxP,buf,buf_ub);
    NTL_Pack(x.auxQ,buf,buf_ub);
    }
  };

bool EC_Unpack(ec_point& x, unsigned char * buf, int& buf_ub, int buf_size) {
  if (buf_ub>=buf_size) return false;
  if (buf[buf_ub]==0) {
    buf_ub++;
    x=to_ec_point(true);
    }
  else {
    buf_ub++;
    ZZ ec_x,ec_y,ec_auxP,ec_auxQ;
    NTL_Unpack(ec_x,buf,buf_ub,buf_size);
    NTL_Unpack(ec_y,buf,buf_ub,buf_size);
    NTL_Unpack(ec_auxP,buf,buf_ub,buf_size);
    NTL_Unpack(ec_auxQ,buf,buf_ub,buf_size);
    x=to_ec_point(ec_x,ec_y);
    x.auxP=ec_auxP;
    x.auxQ=ec_auxQ;
    }
  return true;
  };
