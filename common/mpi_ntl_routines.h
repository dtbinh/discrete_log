#ifndef MPI_NTL_ROUTINES_H
#define MPI_NTL_ROUTINES_H 1

#include "mpi.h"
#include <NTL/ZZ.h>

NTL_CLIENT

#define MPI_NTL_SIGN 5000
#define MPI_NTL_DATA 5001

#define MPI_BUF_LEN  5100
#define MPI_BUF_DATA 5101

//Структура для передачи параметров NTL-числа
struct ntl_params {
int sign; //Знак
int size; //Длина двоичного представления числа в байтах. 0 - число равно 0.
};

//Послать NTL число
int MPI_NTL_Send(ZZ  x, int process_id);
//Принять NTL число
int MPI_NTL_Recv(ZZ& x, int process_id);
//Широковещательно разослать NTL число
int MPI_NTL_Bcast(ZZ& x, int process_id);

//Функции для групповой работы
void NTL_Pack(ZZ x, unsigned char * buf, int& buf_ub);
bool NTL_Unpack(ZZ& x, unsigned char * buf, int& buf_ub, int buf_size);

void Buf_send(unsigned char* buf, int buf_size, int process_id);
void Buf_recv(unsigned char* buf, int& buf_size, int process_id);
void Buf_bcast(unsigned char* buf, int& buf_size, int process_id);

#endif
