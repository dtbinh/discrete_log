#ifndef EC_POINT_H
#define EC_POINT_H 1

#include <NTL/ZZ.h>
#include <iostream>
NTL_CLIENT

#pragma once

int ec_compare(const void * ec1, const void * ec2);
int ec_compare_bs(const void * ec1, const void * ec2);

//Класс, хранящий параметры нашей эллиптической кривой
class ec_curve {
		//Эллиптическая кривая задаётся уравнением
		//y^2==x^3+A*x+B
		//Если х-ка поля !=2 и !=3
	public: ZZ A,B;
		//Модуль
		ZZ p;
		static ec_curve * Instance();
	private:
//		static ec_curve * _instance;
	};

//ec_curve * ec_curve::_instance=NULL;

class ec_point {

public:

//Вспомогательные переменные, используются при хранении скалярных множителей P,Q
ZZ auxP,auxQ;

//Свойства
private:

//Координаты, однозначно определяющие точку
ZZ x,y;

//Является ли точка бесконечно удалённой (нулём аддитивной группы)
bool is_inf;

public:
//Конструктор
ec_point();
ec_point(ZZ x, ZZ y);
ec_point(bool is_inf);
ec_point(const ec_point& op2);

static void Init(ZZ A, ZZ B, ZZ p);
static ZZ get_p();
//Вывод на печать
void print(ostream& out_stream);

//Является ли бесконечностью
bool is_infinity();

//Обратный эл-т по сложению
ec_point operator -();

//Сумма 2-х точек
ec_point operator +(ec_point op2);

bool operator ==(ec_point op2);
//Умножение точки на число
ec_point double_point();
ec_point operator *(ZZ op2);

//Отношения порядка
int friend ec_compare(ec_point ec1, ec_point ec2);
int friend ec_compare_bs(ec_point ec1, ec_point ec2);

//Получение координат
ZZ get_x();
ZZ get_y();
};

ec_point to_ec_point(ZZ x, ZZ y);
ec_point to_ec_point(bool is_inf);
#endif
