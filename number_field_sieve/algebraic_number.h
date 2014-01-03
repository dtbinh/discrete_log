#ifndef ALGEBRAIC_NUMBER
#define ALBEBRAIC_NUMBER 1
#include <NTL/ZZ.h>

NTL_CLIENT

//Класс для хранения параметров расширенного числового поля
class number_field {
public: ZZ sqr_alpha; //Квадрат alpha
};

//Класс, представляющий собой алгебраическое число
class algebraic_number {
private:
ZZ c,d;  //c+d*alpha

public:
 //Статическая функция, инициализирующая числовое поле
static void init_sqr_alpha(ZZ new_alpha);
//Конструкторы алгебраического числа
algebraic_number(ZZ c,ZZ d);
algebraic_number(long c,long d);
algebraic_number(); //0+0*alpha

//Получение нормы числа
ZZ norm();

//Арифметические операции
algebraic_number operator + (algebraic_number op2);
algebraic_number operator - (algebraic_number op2);
algebraic_number operator * (algebraic_number op2);
algebraic_number operator * (ZZ op2);
algebraic_number operator / (algebraic_number op2);

//Возведение в заданную степень
algebraic_number power(ZZ x);
algebraic_number power(long x);

//Является ли алгебраическое число делителем числа
bool is_divisible (algebraic_number cand);

//Получить комплексно-сопряжённое число
algebraic_number get_mate();

//Получить c
ZZ get_c();

//Получить d
ZZ get_d();

//Вывести на стандартный вывод
void print();
};

//Переводим в алгебраическое число
algebraic_number to_algebraic_number(ZZ c, ZZ d);

#endif
