#include "pollard_factor.h"
#include <NTL/ZZ.h>
#include <fstream>
using namespace std;
NTL_CLIENT

int main() {
ZZ n;
ZZ x0,q;
long rez_count,i;
ZZ results[MAX_FACT];
ZZ subgroups[MAX_FACT];
long subgroups_count;
ifstream fin("input.txt");
ofstream fout("output.txt");
fin >> n;

pollard_factor(n,results,rez_count);

//find_subgroup_orders(results,rez_count,subgroups,subgroups_count);

for (i=0;i<rez_count;i++) {
	fout << results[i] << endl;
	}
fin.close();
fout.close();
}
