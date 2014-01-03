#include "gauss.h"
#include <NTL/ZZ.h>
NTL_CLIENT

//Обработка строки матрицы по опорному вектору
//Обрабатывается переменная start, считается, что до неё строки занулены
int make_vector(ZZ vec_to_correct[], ZZ base_vec[], int start, int size, ZZ r) {
int i;
ZZ x,y,mul;
mul = vec_to_correct[start];
for (i=start;i<size;i++) {
	x = MulMod(vec_to_correct[i],base_vec[start],r);
	y = MulMod(base_vec[i],mul,r);
	vec_to_correct[i]=SubMod(x,y,r);
	}
};

//Обработка строки матрицы по опорному вектору
//Никаких предположений по поводу зануления строк не делается
int make_common_vector(ZZ vec_to_correct[], ZZ base_vec[], int start, int size, ZZ r) {
int i;
ZZ x,y,mul;
mul = vec_to_correct[start];

if (r==0) {
	for (i=0;i<size;i++) {
		x = vec_to_correct[i]*base_vec[start];
		y = base_vec[i]*mul;
		vec_to_correct[i]=x-y;
		}
	}
else {
	for (i=0;i<size;i++) {
		x = MulMod(vec_to_correct[i],base_vec[start],r);
		y = MulMod(base_vec[i],mul,r);
		vec_to_correct[i]=SubMod(x,y,r);
		}
	}
};
