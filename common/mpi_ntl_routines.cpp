#include "mpi_ntl_routines.h"
#include "mpi.h"
#include <NTL/ZZ.h>

NTL_CLIENT

//Функция отсылки ZZ числа
int MPI_NTL_Send(ZZ  x, int process_id) {
unsigned char * buf;
ntl_params def_data;
def_data.size = NumBytes(x);
def_data.sign = x>=0?1:-1;
MPI_Send(&def_data,sizeof(ntl_params),MPI_CHAR,process_id,MPI_NTL_SIGN,MPI_COMM_WORLD);

if (def_data.size > 0) {
	buf = new unsigned char [def_data.size];
	BytesFromZZ(buf,x,def_data.size);
	MPI_Send(buf,def_data.size,MPI_CHAR,process_id,MPI_NTL_DATA,MPI_COMM_WORLD);
	delete [] buf;
	}
return 0;
};

//Функция приёма ZZ числа
int MPI_NTL_Recv(ZZ& x, int process_id) {
MPI_Status status;
unsigned char * buf;
ntl_params def_data;
MPI_Recv(&def_data,sizeof(ntl_params),MPI_CHAR,process_id,MPI_NTL_SIGN,MPI_COMM_WORLD,&status);
if (def_data.size>0) {
	buf = new unsigned char [def_data.size];
	MPI_Recv(buf,def_data.size,MPI_CHAR,process_id,MPI_NTL_DATA,MPI_COMM_WORLD,&status);
	ZZFromBytes(x,buf,def_data.size);
	x*=def_data.sign;
	delete [] buf;
	}
else {
	x=0;
	}
return 0;
};

//Функция широковещательной рассылки ZZ числа
int MPI_NTL_Bcast(ZZ& x, int process_id){
unsigned char * buf;
ntl_params def_data;
def_data.size = NumBytes(x);
def_data.sign = x>=0?1:-1;

MPI_Bcast(&def_data,sizeof(ntl_params),MPI_CHAR,process_id,MPI_COMM_WORLD);

if  (def_data.size >0) {
	buf = new unsigned char [def_data.size];
	BytesFromZZ(buf,x,def_data.size);
	MPI_Bcast(buf,def_data.size,MPI_CHAR,process_id,MPI_COMM_WORLD);
	ZZFromBytes(x,buf,def_data.size);
	x*=def_data.sign;
	delete [] buf;
	}
else {
	x=0;
	}
return 0;
};


//Фунции для групповой работы
void NTL_Pack(ZZ x, unsigned char * buf, int& buf_ub) {
  ntl_params def_data;
  def_data.size = NumBytes(x);
  def_data.sign = x>=0?1:-1;
  memcpy(&(buf[buf_ub]),&def_data,sizeof(def_data));
  buf_ub+=sizeof(def_data);
  if (def_data.size>0) {
    unsigned char * temp_buf = new unsigned char [def_data.size];
    BytesFromZZ(temp_buf,x,def_data.size);  
    memcpy(&(buf[buf_ub]),temp_buf,def_data.size);
    buf_ub+=def_data.size;
    delete [] temp_buf;
    }
  };

bool NTL_Unpack(ZZ& x, unsigned char * buf, int& buf_ub, int buf_size){
  if (buf_ub>=buf_size) return false;
  ntl_params def_data;
  memcpy(&def_data,&(buf[buf_ub]),sizeof(def_data));
  buf_ub+=sizeof(def_data);
  if (def_data.size>0) {
    unsigned char * temp_buf = new unsigned char [def_data.size];
    memcpy(temp_buf,&(buf[buf_ub]),def_data.size);
    ZZFromBytes(x,temp_buf,def_data.size);
    buf_ub+=def_data.size;
    delete [] temp_buf;
    }
  else {
    x=0;
    }
  x*=def_data.sign;
  return true;
  };

void Buf_send(unsigned char * buf, int buf_size, int process_id) {
  MPI_Send(&buf_size,1,MPI_INT,process_id,MPI_BUF_LEN,MPI_COMM_WORLD);
  MPI_Send(buf,buf_size,MPI_CHAR,process_id,MPI_BUF_DATA,MPI_COMM_WORLD);
  }

void Buf_recv(unsigned char * buf, int& buf_size, int process_id) {
  MPI_Status status;
  MPI_Recv(&buf_size,1,MPI_INT,process_id,MPI_BUF_LEN,MPI_COMM_WORLD,&status);
  MPI_Recv(buf,buf_size,MPI_CHAR,process_id,MPI_BUF_DATA,MPI_COMM_WORLD,&status);
  }

void Buf_bcast(unsigned char * buf, int& buf_size, int process_id) {
  MPI_Bcast(&buf_size,1,MPI_INT,process_id,MPI_COMM_WORLD);
  MPI_Bcast(buf,buf_size,MPI_CHAR,process_id,MPI_COMM_WORLD);
  }

#ifdef MODULE_TEST

int main() {
ZZ x=to_ZZ(5);
ZZ y=to_ZZ(0);
ZZ z=to_ZZ(-3);

unsigned char data[1000];
int data_len;
int data_ub;

data_len=0;
cout << "x=" << x << endl;
cout << "y=" << y << endl;
cout << "z=" << z << endl;

NTL_Pack(x,data,data_len);
NTL_Pack(y,data,data_len);
NTL_Pack(z,data,data_len);
cout << "data_len = " << data_len << endl;

data_ub=0;
NTL_Unpack(x,data,data_ub,data_len);
NTL_Unpack(y,data,data_ub,data_len);
NTL_Unpack(z,data,data_ub,data_len);
cout << "x=" << x << endl;
cout << "y=" << y << endl;
cout << "z=" << z << endl;

cout << "data_ub = " << data_ub << endl;
}

#endif
