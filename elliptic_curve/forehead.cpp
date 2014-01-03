#include "ec_point.h"
#include "ec_routine.h"
#include <time.h>
#include <NTL/ZZ.h>
#include <iostream>
#include <fstream>
using namespace std;


int main() {
ec_point P,Q;
ZZ r;
//Чтение исходных данных
//Ввести данные
//Кривую
ifstream ftask("task.txt");
ZZ A,B;
ftask >> A >> B;
//Модуль
ZZ p;
ftask >> p;
ec_point::Init(A,B,p);

//Образующую точку
ZZ x,y;
ftask >> x >> y;
Q= to_ec_point(x,y);
//Порядок точки
ftask >> r;
//Сама точка
ftask >> x >> y;
P=to_ec_point(x,y);
ftask.close();

ZZ i;
P=Q;
ZZ c=to_ZZ(1);
for (i=to_ZZ(0);;i++) {
	if (P.is_infinity()) break;
	c++;
	P=P+Q;
	}
cout << c << endl;
return 0;
}
