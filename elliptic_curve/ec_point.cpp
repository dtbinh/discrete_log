#include "ec_point.h"
#include <NTL/ZZ.h>
#include <iostream>
using namespace std;
NTL_CLIENT

ec_point to_ec_point(ZZ x, ZZ y) {
ec_point Q(x,y);
return Q;
};

ec_point to_ec_point(bool is_inf) {
ec_point Q(is_inf);
return Q;
};

//Данная реализация синглетона подходит только для однопоточных программ
ec_curve * ec_curve::Instance() {
static ec_curve inst;
return &inst;
};

void ec_point::Init(ZZ _A, ZZ _B, ZZ _p) {
ec_curve * ec_singleton = ec_curve::Instance();
ec_singleton->A = _A;
ec_singleton->B = _B;
ec_singleton->p = _p;
}

//Конструктор
ec_point::ec_point(){
x = to_ZZ(0);
y = to_ZZ(0);
is_inf = true;
};

ec_point::ec_point(ZZ _x, ZZ _y) {
x = _x;
y = _y;
is_inf = false;
};

ec_point::ec_point(bool _is_inf) {
is_inf = _is_inf;
}

ec_point::ec_point(const ec_point& op2){
x = op2.x;
y = op2.y;
auxP = op2.auxP;
auxQ = op2.auxQ;
is_inf = op2.is_inf;
};

//Вывод на печать
void ec_point::print(ostream& out_stream) {
if (is_inf)
	out_stream << "INF" << endl;
else 
	out_stream << x << " " << y << endl;
}

//Является ли бесконечностью
bool ec_point::is_infinity(){
return is_inf;
};

//Обратный эл-т по сложению
ec_point ec_point::operator -(){
ec_curve * ec_singleton = ec_curve::Instance();
if (is_inf)
	return ec_point(true);
else 
	return ec_point(x,NegateMod(y,ec_singleton->p));
};

//Сумма 2-х точек
//a4=A
//a6=B
//a1,a2,a3=0
ec_point ec_point::operator +(ec_point op2){
	ec_curve * ec_singleton = ec_curve::Instance();
	ZZ x3,y3,lambda,nu;
	//Обработка случаев, если среди операндов нуль группы
	if (is_infinity() && op2.is_infinity()) return ec_point(true);
	if (is_infinity()) return ec_point(op2.x,op2.y);
	if (op2.is_infinity()) return ec_point(x,y);

	//Нулей нет, дальше всё по Вейерштрассу
	if ((x==op2.x) && (y==NegateMod(op2.y,ec_singleton->p))) {
		//Если две противоположные точки, результатом будет нуль группы
		return ec_point(true);
		}
	else {
		if (x!=op2.x) { //Если сладываем две разные точки
			lambda = (((op2.y-y)%ec_singleton->p) * InvMod((op2.x-x)%ec_singleton->p,ec_singleton->p))%ec_singleton->p;
			}
		else { //Если удваиваем точку
			lambda = ((3*x*x+ec_singleton->A)%ec_singleton->p * InvMod((2*y)%ec_singleton->p,ec_singleton->p))%ec_singleton->p;
			}
		x3 = (lambda*lambda - x - op2.x)%ec_singleton->p;
		y3 = (lambda*(x-x3)-y)%ec_singleton->p;
		return ec_point(x3,y3);
		}
	}

//Удвоение точки - можно съэкономить время, избежав лишних проверок
ec_point ec_point::double_point(){
ec_curve * ec_singleton = ec_curve::Instance();
ZZ lambda,x3,y3;
if (is_infinity())
	return ec_point(true);
else {
	lambda = ((3*x*x+ec_singleton->A)%ec_singleton->p * InvMod((2*y)%ec_singleton->p,ec_singleton->p))%ec_singleton->p;
	x3 = (lambda*lambda - to_ZZ(2)*x)%ec_singleton->p;
	y3 = (lambda*(x-x3)-y)%ec_singleton->p;
	return ec_point(x3,y3);
	}
}

//Умножение точки на число
ec_point ec_point::operator *(ZZ op2) {
	ec_point acc(true);
	ec_point mul(x,y);
	int sign;
	if (op2>0) {
		sign=1;
		}
	else if (op2<0) {
		sign=-1;
		}
	else {
		return ec_point(true);
		}

	op2=op2*sign;

	while (op2!=to_ZZ(0)) {
		if (IsOdd(op2))
			acc=acc+mul;
		op2>>=1;
		mul=mul.double_point();
		}
	if (sign==1) {
		return acc;
		}
	else {
		return -acc;
		}
	}

//Порядок задаётся сравнением x-координаты точек
//Если x-координата совпадает, сравниваются y-координаты
//Нуль группы считается меньше любой точки
int ec_compare(ec_point ec1, ec_point ec2) {
	if (!ec1.is_infinity() && !ec2.is_infinity()) {
		if (ec1.x>ec2.x) return 1;
		else if (ec1.x<ec2.x) return -1;
		else if (ec1.x==ec2.x) {
			if (ec1.y>ec2.y) return 1;
			else if (ec1.y<ec2.y) return -1;
			else if (ec1.y==ec2.y) return 0;
			}
		}
	else if (ec1.is_infinity() && ec2.is_infinity())
		return 0;
	else if (ec1.is_infinity() && !ec2.is_infinity())
		return -1;
	else if (!ec1.is_infinity() && ec2.is_infinity())
		return 1;

	}

//Аналог для использования в функции qsort
int ec_compare(const void * ec1, const void * ec2) {
	return ec_compare(*(ec_point*)ec1,*(ec_point*)ec2);
	}

//Функция для использования в бинарном поиске
//Считает одинаковыми точку и обратную ей ( x-координаты совпадают)
//Порядок задаётся сравнением x-координаты точек
//Нуль группы считается меньше любой точки
int ec_compare_bs(ec_point ec1, ec_point ec2) {
	if (!ec1.is_infinity() && !ec2.is_infinity()) {
		if (ec1.x>ec2.x) return 1;
		else if (ec1.x<ec2.x) return -1;
		else if (ec1.x==ec2.x) {
			return 0;
			}
		}
	else if (ec1.is_infinity() && ec2.is_infinity())
		return 0;
	else if (ec1.is_infinity() && !ec2.is_infinity())
		return -1;
	else if (!ec1.is_infinity() && ec2.is_infinity())
		return 1;

	}

//Аналог для использования в функции qsort
int ec_compare_bs(const void * ec1, const void * ec2) {
	return ec_compare_bs(*(ec_point*)ec1,*(ec_point*)ec2);
	}


//Получение координат
ZZ ec_point::get_x(){
	return x;
	};
ZZ ec_point::get_y(){
	return y;
	};

ZZ ec_point::get_p() {
	ec_curve * ec_singleton = ec_curve::Instance();
	return ec_singleton->p;
	}

bool ec_point::operator ==(ec_point op2){
if (x==op2.x && y==op2.y)
	return true;
else
	return false;
};
