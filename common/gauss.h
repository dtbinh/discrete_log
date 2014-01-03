#ifndef GAUSS_H
#define GAUSS_H 1

#include <NTL/ZZ.h>
NTL_CLIENT

int make_vector(ZZ vec_to_correct[], ZZ base_vec[], int start, int size, ZZ r);
//Обработка строки матрицы по опорному вектору
//Никаких предположений по поводу зануления строк не делается
int make_common_vector(ZZ vec_to_correct[], ZZ base_vec[], int start, int size, ZZ r);
#endif
