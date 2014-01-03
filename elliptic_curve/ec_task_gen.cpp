#include "ec_point.h"
#include "ec_routine.h"
#include "pollard_factor.h"
#include <time.h>

#include <NTL/ZZ.h>
#include <iostream>
using namespace std;

int main() {
ZZ p;
ZZ A,B;
long p_len;

SetSeed(to_ZZ(time(NULL)));
//Получить длину модуля
cin >> p_len;
//Получить модуль как простое число заданной длины
p = GenPrime_ZZ(p_len);

//Сформировать кривую
A=to_ZZ(1);
B=to_ZZ(1);

ec_point::Init(A,B,p);

//Найти образующую точку
ec_point Q;
ZZ x,y;

x=RandomBnd(p/2);
while (Jacobi((x*x*x+A*x+B)%p,p)!=1)
	x++;
y=SqrRootMod((x*x*x+A*x+B)%p,p);

Q = to_ec_point(x,y);

//Найти порядок точки
ZZ r = Shanks_cycle_group_order(Q);


/*
//Получить возможные порядки подгрупп
ZZ tmp_vec[MAX_FACT],rez_vec[MAX_FACT];
long tmp_ub, rez_ub;

pollard_factor(r, tmp_vec, tmp_ub);
find_subgroup_orders(tmp_vec, tmp_ub, rez_vec, rez_ub);

bool flag;
long i;

cerr << "==============================" << endl;
for (i=0;i<rez_ub;i++)
	cerr << rez_vec[i] << " ";
cerr << endl;
cerr << "==============================" << endl;
//do {
	
	//flag = false;
	for (i=0;i<rez_ub;i++)
		if ((Q*rez_vec[i]).is_infinity()) {
			r=rez_vec[i];
			}
//	} while (flag);
*/

ZZ l = RandomBnd(r-2)+2;
ec_point P = Q*l;

//Вывести данные
//Кривую
cout << A << " " << B << endl;
//Модуль
cout << p << endl;
//Образующую точку
Q.print(cout);
//Порядок точки
cout << r << endl;
P.print(cout);
return 0;
}
