#include "ec_routine.h"

//Осуществим бинарный поиск в БД по x-координате точки
long find_by_x(ec_point * data, long lb, long rb, ZZ x) {
int rez;
while (lb<=rb) {
	rez = (lb+rb)>>1;
	if (data[rez].get_x()==x) {
		return rez;
		}
	else if (data[rez].get_x()<x) {
		lb=rez+1;
		}
	else {
		rb=rez-1;
		}
	}
return -1;
}

//Реализуем алгоритм Шенкса для нахождения порядка циклической группы
ZZ Shanks_cycle_group_order(ec_point Q) {

ZZ p = ec_point::get_p();

ec_point * db = NULL;
long k = trunc_long(SqrRoot(SqrRoot(2*p)),32)+2;
long i;
db = new ec_point[k];
db[0]=Q;
db[0].auxQ=to_ZZ(1);


for (i=1;i<k;i++) {
	db[i]=db[i-1]+Q;
	db[i].auxQ=db[i-1].auxQ+1;
	}

qsort(db,k,sizeof(ec_point),ec_compare);

ec_point P,R;
P = Q*to_ZZ(2*k+1);
R = Q*to_ZZ(p+1);

ec_point acc(true),test;
ZZ d,e,t;
long ind=-1;

for (i=0;i<k;i++) {
	test=R+acc;
	ind=find_by_x(db,0,k-1,test.get_x());
	if (ind!=-1) {
		t=i;
		break;
		}
	test=R+-acc;
	ind=find_by_x(db,0,k-1,test.get_x());
	if (ind!=-1){
		t=-i;
		break;
		}
	acc=acc+P;
	}

ZZ n;
if (ind==-1) {
	return to_ZZ(-1);
	}
else {
	e=db[ind].auxQ;
	if (ec_compare(test,db[ind])) e*=-1;
	d=2*t;
	e=e-t;
	n = p+1+d*k-e;
	return n;
	}
}

void init_round_robin_scheme(long csize,long**& round_robin_scheme, long& rounds) {
long * half1 = new long [csize+5];
long * half2 = new long [csize+5];
int i,j,t;
long tsize = csize;
if (csize%2==1) tsize++;
for (i=0;i<tsize/2;i++) {
  half1[i]=i;
  half2[i]=i+tsize/2;
  }
if (csize!=tsize)
  half2[tsize/2-1]=-1;
round_robin_scheme = new long*[csize];
rounds=tsize-1;
for (i=0;i<csize;i++)
  round_robin_scheme[i] = new long[rounds];
for (i=0;i<rounds;i++) {
  //Заполним матрицу round_robin_scheme
  for (j=0;j<tsize/2;j++) {
    if (half1[j]!=-1)
      round_robin_scheme[half1[j]][i]=half2[j];
    if (half2[j]!=-1)
      round_robin_scheme[half2[j]][i]=half1[j];
    }
  
  //Произведём сдвиг
  
  t=half2[0];
  for (j=0;j<(tsize/2)-1;j++)
    half2[j]=half2[j+1];
  half2[(tsize/2)-1]=half1[(tsize/2)-1];
  
  for (j=(tsize/2)-1;j>1;j--)
    half1[j]=half1[j-1];
  half1[1]=t;
  
  }  

} 

ZZ calc_answer(ec_point p1, ec_point p2, ec_point Q, ec_point P, ZZ r) {

ZZ b=p1.auxP;
ZZ a=p1.auxQ;

ZZ d=p2.auxP;
ZZ c=p2.auxQ;

if (p1.get_y()!=p2.get_y()) {
	//if (c!=0)
	  c=NegateMod(c,r);
	//if (d!=0)
	  d=NegateMod(d,r);
	}

ZZ l=(c-a)%r;
ZZ inv_part=(b-d)%r;

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

return l;
}
