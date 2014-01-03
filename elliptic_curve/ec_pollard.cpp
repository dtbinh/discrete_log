#include "ec_point.h"
#include "ec_routine.h"
#include <time.h>
#include <NTL/ZZ.h>
#include <iostream>
#include <fstream>
#include "mpi.h"
using namespace std;

ZZ one_third, two_third;
ec_point P,Q;
ZZ r;

void tau(ec_point& X) {
ZZ auxQ,auxP;
if (X.get_x()<one_third) {
	auxQ=X.auxQ+1;
	auxP=X.auxP;
	X=X+Q;
	X.auxQ = auxQ%r;
	X.auxP = auxP%r;
	}
else if (X.get_x()>one_third && X.get_x()<two_third) {
	auxQ=X.auxQ;
	auxP=X.auxP+1;
	X=X+P;
	X.auxQ = auxQ%r;
	X.auxP = auxP%r;
	}
else {
	auxQ=X.auxQ+2;
	auxP=X.auxP;
	X=X+Q*to_ZZ(2);
	X.auxQ = auxQ%r;
	X.auxP = auxP%r;
	}
}

int main(int argc, char ** argv) {
MPI_Init(&argc,&argv);

SetSeed(to_ZZ(time(NULL)));
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
//Задача P=Q*l
//Найти l (mod r)

//Устанавливаем параметры для отображения tau
one_third = p/3;
two_third = one_third*2;

//Инициализируем R,S
ZZ a=RandomBnd(r),b=RandomBnd(r),c=RandomBnd(r),d=RandomBnd(r);

/*
ec_point R=(Q*a) + (P*b);
R.aux1=a;
R.aux2=b;

ec_point S=(Q*c) + (P*d);
S.aux1=c;
S.aux2=d;
*/

ec_point R=Q;
R.auxQ=1;
R.auxP=0;

ec_point S=Q;
S.auxQ=1;
S.auxP=0;


double t1,t2;
t1=MPI_Wtime();

//Начинаем цикл метода Полларда
do {
	tau(R);
	tau(S);
	tau(S);
	} while (R.get_x()!=S.get_x());

t2=MPI_Wtime();
cout << "Time: " << t2-t1 << endl;

a=R.auxQ;
b=R.auxP;
c=S.auxQ;
d=S.auxP;

if (R.get_y()!=S.get_y()) {
	c=NegateMod(c,r);
	d=NegateMod(d,r);
	}

ZZ l=(c-a)%r;
ZZ inv_part=b-d;

	ZZ divisor = GCD(inv_part,r);

	ZZ tr=r/divisor;
	ZZ t_norm;
	//Правило гласит, что a=b (mod rs), тогда и только тогда, когда a=b (mod r) и a=b (mod s)
	//При этом НОД(r,s)=1, иначе правило работает не всегда
	do {
		t_norm=GCD(divisor,tr);
		tr/=t_norm;
		divisor*=t_norm;
		} while (t_norm!=1);

	cout << "divisor of r = " << divisor << endl;
	l = MulMod(InvMod(inv_part%tr,tr),l%tr,tr);

	//А теперь скорректируем ответ, если divisor!=1
	ZZ corr;
	for (corr=0;corr<r;corr+=tr) {
		if (Q*(l+corr)==P) {
			l+=corr;
			break;
			}
		}
	l%=r;
	cout << "Answer = " << l << endl;
	//Проверим, правильный ли мы получили ответ
	if (Q*l==P) {
		cout << "Answer is correct" << endl;
		}
	else {
		cout << "Answer is WRONG!" << endl;
		}

MPI_Finalize();
return 0;
}
