#include "algebraic_number.h"
#include <iostream>
using namespace std;

//Класс для хранения параметров числового поля
number_field nf_singleton;

//Инициализация числового поля
void algebraic_number::init_sqr_alpha(ZZ new_alpha){
nf_singleton.sqr_alpha = new_alpha;
}

//Конструкторы
algebraic_number::algebraic_number(ZZ _c,ZZ _d) {
c=_c;
d=_d;
};

algebraic_number::algebraic_number(long _c,long _d){
c=to_ZZ(_c);
d=to_ZZ(_d);
};

algebraic_number::algebraic_number() {
c=0;
d=0;
};

//Получение нормы алгебраического числа
ZZ algebraic_number::norm() {
return c*c-d*d*nf_singleton.sqr_alpha;
};

//Реализация арифметических операций
algebraic_number algebraic_number::operator + (algebraic_number op2){
algebraic_number rez;
rez.c = c+op2.c;
rez.d = d+op2.d;
return rez;
};

algebraic_number algebraic_number::operator - (algebraic_number op2){
algebraic_number rez;
rez.c = c-op2.c;
rez.d = d-op2.d;
return rez;
};

algebraic_number algebraic_number::operator * (algebraic_number op2){
algebraic_number rez;
rez.c = c*op2.c + d*op2.d*nf_singleton.sqr_alpha;
rez.d = c*op2.d + d*op2.c;
return rez;
};

algebraic_number algebraic_number::operator * (ZZ op2) {
algebraic_number rez;
rez.c = c*op2;
rez.d = d*op2;
return rez;
};

//Быстрое возведение алгебраического числа в степень
algebraic_number algebraic_number::power(ZZ x){
algebraic_number rez(1,0);
algebraic_number acc=*this;
while (x>0) {
	if (x%2==1) {
		rez=rez*acc;
		}
	x/=2;
	acc=acc*acc;
	}
return rez;
};

algebraic_number algebraic_number::power(long x) {
return power(to_ZZ(x));
};

//Деление реализовано домножением числителя и знаменателя на комплексно-сопряжённое число
algebraic_number algebraic_number::operator / (algebraic_number op2){
algebraic_number rez;
ZZ n;
rez.c = c*op2.c - d*op2.d*nf_singleton.sqr_alpha;
rez.d = d*op2.c - c*op2.d;
n=op2.norm();
rez.c /=n;
rez.d /=n;
return rez;
};


//Делится ли число на число
bool algebraic_number::is_divisible (algebraic_number cand){
if (norm()%cand.norm()!=0)
	return false;
algebraic_number rez;
ZZ n;
rez.c = c*cand.c - d*cand.d*nf_singleton.sqr_alpha;
rez.d = c*cand.d - d*cand.c;
n=cand.norm();
if ((rez.c % n == 0) && (rez.d % n == 0)) {
	return true;
	}
else {
	return false;
	}
};

//Вывод на стандартный вывод
void algebraic_number::print(){
cout << c << " " << d << endl;
};

//перевод в алгебраическое число
algebraic_number to_algebraic_number(ZZ c, ZZ d){
algebraic_number rez(c,d);
return rez;
};

//Получить комплексно-сопряжённое число
algebraic_number algebraic_number::get_mate(){
algebraic_number rez(c,-d);
return rez;
};

//Получить c
ZZ algebraic_number::get_c(){
return c;
};

//Получить d
ZZ algebraic_number::get_d(){
return d;
};
