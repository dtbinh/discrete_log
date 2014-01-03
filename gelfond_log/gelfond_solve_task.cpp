#include <NTL/ZZ.h>
#include <iostream>
#include <fstream>
using namespace std;

NTL_CLIENT
#define MAX_UB 100
//Функция поиска логарифма
//Поскольку r невелико, используем полный перебор
ZZ find_log(ZZ a, ZZ b, ZZ p,ZZ r) {
ZZ i;
for (i=to_ZZ(0);i<r;i++) //Пробуем все возможные логарифмы
	if (PowerMod(a,i,p)==b) {
		return i; //Нашли результат
		}
return to_ZZ(-1); //Ошибка во входных данных
}

int main() {
ifstream ftask("task.txt");
ZZ a,b,p,r; //Исходные данные задачи
long i,j,k;
ftask >> a >> b >> p >> r;

ZZ prime_vec[MAX_UB]; //Вектор простых сомножителей
long index_vec[MAX_UB]; //Вектор их показателей
long prime_ub; //Число простых сомножителей
ftask >> prime_ub;
for (i=0;i<prime_ub;i++)
  ftask >> prime_vec[i] >> index_vec[i];
ftask.close();

//Считали данные
//Теперь решаем задачу
//Находим частичные логарифмы
ZZ l[MAX_UB];
for (i=0;i<prime_ub;i++) {
  ZZ x = prime_vec[i];
  l[i]=0;
  ZZ ta,tb=to_ZZ(0),S=to_ZZ(1);
  ta = PowerMod(a,r/x,p);

  for (j=0;j<index_vec[i];j++) {
    tb=PowerMod(MulMod(b,InvMod(S,p),p),r/x,p);
    ZZ k=find_log(ta,tb,p,x);
    cout << ta << " " << tb << " " << x << " " << k << endl;
    l[i]+=k*(x/prime_vec[i]);
    S = (S*PowerMod(a,k*(x/prime_vec[i]),p))%p;
    x*=prime_vec[i];
    }
  }

//Восстанавливаем исходный логарифм по китайской теореме об остатках
ZZ ans,ans1,tp,tp1;
ans=l[0];
tp=power(prime_vec[0],index_vec[0]);
for (i=1;i<prime_ub;i++) {
  ans1=l[i];
  tp1=power(prime_vec[i],index_vec[i]);
  CRT(ans1,tp1,ans,tp);
  ans=(ans1+tp1)%tp1; //Функция может вернуть отрицательное значение, что не есть хорошо
  tp=tp1;
  }

//Вывод ответа
cout << ans << endl;

//Проверка результата
if (PowerMod(a,ans,p)==b) {
  cout << "Right" << endl;
  }
else {
  cout << "WRONG" << endl;
  }
return 0;
}
