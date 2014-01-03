#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

NTL_CLIENT

int main() {
long x;
cin >> x;
cout << RandomPrime_ZZ(x)*RandomPrime_ZZ(x) << endl;
}
