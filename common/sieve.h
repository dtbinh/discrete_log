#ifndef SIEVE_H
#define SIEVE_H 1
#include <NTL/ZZ.h>
NTL_CLIENT

//Общая функция просеивания
bool sieving(ZZ cand, ZZ p, ZZ basis_mul, ZZ basis[],int basis_size, ZZ result_powers[]);
//Наивное просеивание - пытаемся делить на элементы базиса
bool naive_sieving(ZZ cand, ZZ p, ZZ basis[],int basis_size, ZZ result_powers[]);
//Просеивание поиском НОД числа и произведения эл-тов базиса, взятых с определёнными степенями
bool gcd_sieving(ZZ cand, ZZ p, ZZ basis_mul, ZZ basis[],int basis_size, ZZ result_powers[]);
#endif
