#include <NTL/ZZ.h>
#include <fstream>
NTL_CLIENT

int main() {
//основание, степень, модуль, порядок группы
ZZ a,b,p,r; 
ifstream ftask("task.txt");
//Вводим исходные данные
ftask >> a >> b >> p >> r;
ftask.close();
ZZ i;

for (i=to_ZZ(0);i<r;i++)
	if (PowerMod(a,i,p)==b) {
		cout << "x=" << i << endl;
		return 0;
		}
cout << "NO ANSWER" << endl;
return 0;
}
