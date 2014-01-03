#include "sieve.h"
#define SIEVE_TYPE GCD
//Данный файл содержит функции определения B-гладкости целого числа
//(функции решета)

/**
Абстрактная функция - вызывает одну из двух реализаций
@param cand - Кандидат на гладкость
@param p - модуль, по которому ведутся вычисления
@param basis_mul - Заранее вычисленное число для нахождения НОД
@param basis - Базис разложения
@param basis_size - размер базиса
@param result_powers - степени разложения (если кандидат был разложен)
@return истина - кандидат гладкий, и в result_powers - его разложение. Ложь - кандидат не гладкий
*/

bool sieving(ZZ cand, ZZ p, ZZ basis_mul, ZZ basis[],int basis_size, ZZ result_powers[]) {
#if SIEVE_TYPE == GCD
return gcd_sieving(cand,p,basis_mul,basis,basis_size, result_powers);
#else
return naive_sieving(cand,p,basis,basis_size,result_powers);
#endif
}

/**
В данной функции проверка на гладкость производится поиском НОД с числом,
составленным перемножением чисел базиса, взятых с определённой степенью.
Работает примерно в два раза быстрее проверки "по определению".
@param cand - Кандидат на гладкость
@param p - модуль, по которому ведутся вычисления
@param basis_mul - Заранее вычисленное число для нахождения НОД
@param basis - Базис разложения
@param basis_size - размер базиса
@param result_powers - степени разложения (если кандидат был разложен)
@return истина - кандидат гладкий, и в result_powers - его разложение. Ложь - кандидат не гладкий
*/

bool gcd_sieving(ZZ cand, ZZ p, ZZ basis_mul, ZZ basis[],int basis_size, ZZ result_powers[]) {
int i,cur_div;
ZZ q,r;
//Инициализируем вектор разложения
for (i=0;i<basis_size;i++)
	result_powers[i]=0;

//Если НОД кандидата и проверочного числа равен кандидату,
//то кандидат - гладкий
if (GCD(cand,basis_mul)==cand) {
	cur_div=1;

	//Обычная процедура разложения по базису
	while ((cand>1)&&(cur_div<basis_size)) {
		DivRem(q,r,cand,basis[cur_div]);
		if (r==0) {
			result_powers[cur_div]++;
			cand=q;
			}
		else {
			cur_div++;
			}
	}
	if (cand==1) return true; else return false;
	}
//Если не получилось - стоит попробовать использовать -1 из базиса
cand = NegateMod(cand,p);
result_powers[0]=1;

//Проверяем на гладкость нахождением НОД
if (GCD(cand,basis_mul)==cand) {
	cur_div=1;
	//Гладкий - раскладываем по базису
	while ((cand>1)&&(cur_div<basis_size)) {
		DivRem(q,r,cand,basis[cur_div]);
		if (r==0) {
			result_powers[cur_div]++;
			cand=q;
			}
		else {
			cur_div++;
			}
	}
	if (cand==1) return true; else return false;
	}
return false;
}

/**
В данной функции проверка на гладкость осуществляется "по определению"
Мы просто пытаемся разложить кандидата по базису.
@param cand - Кандидат на гладкость
@param p - модуль, по которому ведутся вычисления
@param basis_mul - Заранее вычисленное число для нахождения НОД
@param basis - Базис разложения
@param basis_size - размер базиса
@param result_powers - степени разложения (если кандидат был разложен)
@return истина - кандидат гладкий, и в result_powers - его разложение. Ложь - кандидат не гладкий
*/
bool naive_sieving(ZZ cand, ZZ p, ZZ basis[],int basis_size, ZZ result_powers[]) {
int i,cur_div=1;
ZZ q,r;

//Инициализируем вектор разложения
for (i=0;i<basis_size;i++)
	result_powers[i]=0;

//Сохраним кандидата на случай
ZZ cand_back=cand;
//Проверяем гладкость разложением по базису
while ((cand>1)&&(cur_div<basis_size)) {
	DivRem(q,r,cand,basis[cur_div]);
	if (r==0) {
		result_powers[cur_div]++;
		cand=q;
		}
	else {
		cur_div++;
		}
	}

//Разложили - отлично
if (cand==1) return true;
//Нет - пробуем применить -1
else {
	cand=cand_back;
	for (i=0;i<basis_size;i++)
		result_powers[i]=0;
	result_powers[0]=1;
	//Умножаем на -1
	cand=NegateMod(cand,p);
	cur_div=1;
	//И опять пробуем разложить
	while ((cand>1)&&(cur_div<basis_size)) {
		DivRem(q,r,cand,basis[cur_div]);
		if (r==0) {
			result_powers[cur_div]++;
			cand=q;
			}
		else {
			cur_div++;
			}
	}
	//Разложили
	if (cand==1) return true;
	//Не разложили
	else return false;
	}
}
