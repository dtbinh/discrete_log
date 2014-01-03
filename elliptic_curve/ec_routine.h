#ifndef EC_ROUTINE_H
#define EC_ROUTINE_H 1
#include "ec_point.h"
#include <iostream>
using namespace std;
//Осуществим бинарный поиск в БД по x-координате точки
long find_by_x(ec_point * data, long lb, long rb, ZZ x);

//Реализуем алгоритм Шенкса для нахождения порядка циклической группы
ZZ Shanks_cycle_group_order(ec_point Q);

//Инициализируем схему обмена информацией
void init_round_robin_scheme(long csize,long**& round_robin_scheme, long& rounds);

//Вычислим ответ по найденной паре
ZZ calc_answer(ec_point p1, ec_point p2, ec_point Q, ec_point P, ZZ r);
#endif
