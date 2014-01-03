#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

NTL_CLIENT

int main() {
ZZ a=to_ZZ(30),b=to_ZZ(10),p=to_ZZ(83);
ZZ x1=to_ZZ("125566390161800"),x2=to_ZZ("674318834107560"),r=to_ZZ(41);
ifstream ftask("task.txt");
ftask >> a >> b >> p >> r;
cout << (PowerMod(a,x1,p)*PowerMod(b,x2,p))%p << endl;
ftask.close();
return 0;
}
