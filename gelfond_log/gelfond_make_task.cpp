#include <iostream>
using namespace std;
#include <NTL/ZZ.h>
#include "pollard_factor.h"
#include <time.h>
NTL_CLIENT

//Составим задачку дискретного логарифмирования
//Которая хорошо решается методом Гельфонда
//Т.е. порядок группы раскладывается на небольшие простые сомножители
int main() {
SetSeed(to_ZZ((int)time(NULL)));
long p_size;
ZZ max_prime;
cin >> p_size >> max_prime;
ZZ p,r;
ZZ results[100];
long rez_count;
long i,j;
bool flag;
do {
  p = RandomPrime_ZZ(p_size);
  r=p-1;
  while (r%2==0) r/=2;
  pollard_factor(r,results,rez_count);
  flag = true;
  for (i=0;i<rez_count;i++)
    if (results[i]>max_prime) 
	flag = false;
  } while (!flag);

r=p-1;

ZZ prime_vec[100];
long index_vec[100];
long prime_ub=0;
for (i=0;i<rez_count;i++) {
    r/=results[i];
    }
prime_vec[prime_ub]=to_ZZ(2);
index_vec[prime_ub]=0;
while (r>1) {
  index_vec[prime_ub]++;
  r/=2;
  }
prime_ub++;

for (i=0;i<rez_count;i++) {
	flag = false;
	for (j=0;j<prime_ub;j++) {
		if (prime_vec[j]==results[i]) {
			flag = true;
			index_vec[j]++;
			break;
			}
		}
	if (!flag) {
		prime_vec[prime_ub]=results[i];
		index_vec[prime_ub]=1;
		prime_ub++;
		}
	}

ZZ a = RandomBnd(p-3)+2;

do {
flag = true;
for (i=0;i<prime_ub;i++) {
  for (j=1;j<=index_vec[i];j++) {
  if (PowerMod(a,(p-1)/power(prime_vec[i],j),p)==1) {
    flag=false;
    a++;
    break;
    }
  }
if (!flag) break;
}
} while (!flag);


ZZ x = RandomBnd(p-3)+2;
cerr << x << endl;
ZZ b = PowerMod(a,x,p);
cout << a << endl;
cout << b << endl;
cout << p << endl;
cout << p-1 << endl;


cout << prime_ub << endl;
for (i=0;i<prime_ub;i++)
  cout << prime_vec[i] << " " << index_vec[i] << endl;
return 0;
}
