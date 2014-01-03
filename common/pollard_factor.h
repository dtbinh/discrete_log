#ifndef POLLARD_FACTOR_H
#define POLLARD_FACTOR_H 1
#include <NTL/ZZ.h>
//Размер вектора факторизуемых элементов
#define MAX_FACT 128
NTL_CLIENT
//Рекурсивная процедура для разложения числа на простые множители алгоритмом Полларда
//в векторе результата будут только простые числа
//их произведение равно исходному числу
//соответственно, одно и то же простое число может встречаться несколько раз
void pollard_factor(ZZ n, ZZ * results, long& rez_count);
//Считая, что факторизованное число было порядком группы
//По разложению на простые числа найдём возможные порядки подгрупп
void find_subgroup_orders(ZZ * input_vec, long input_len, ZZ * output_vec, long & output_len);
#endif
