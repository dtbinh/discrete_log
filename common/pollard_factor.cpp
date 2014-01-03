#include "pollard_factor.h"
//Размер вектора факторизуемых элементов

//Случайное отображение, обладающее сжимающими свойствами
//Для разложения на простые множители алгоритмом Полларда
ZZ f(ZZ x, ZZ n) {
return MulMod(x,x,n);
}

//Рекурсивная процедура для разложения числа на простые множители алгоритмом Полларда
//в векторе результата будут только простые числа
//их произведение равно исходному числу
//соответственно, одно и то же простое число может встречаться несколько раз
void pollard_factor(ZZ n, ZZ * results, long& rez_count) {
	ZZ tmp_rez[MAX_FACT]; //Вектор множителей для найденных сомножителей
	long i,fact_num; //Счётчик цикла, число множителей сомножителя
	ZZ a,b,q; //Параметры алгоритма

	if (ProbPrime(n)) { //Если число простое, это терминальный случай рекурсии
		rez_count=1;
		results[0]=n; //Только один множитель
		return; //Выходим
		}
	else {
		rez_count=0; 
		do {
			a = RandomBnd(n-3)+3; //Инициализируем параметры алгоритма
			b = a;
			do { //Этот цикл, собственно, и реализует метод Полларда
				a = f(a,n);
				b = f(f(b,n),n);
				q = GCD(a-b,n);
			} while (q==1); //Пока не найдём делитель
		} while (q==n); //Найден тривиальный делитель - попробуем ещё разик
				//с другими начальными параметрами
		
		//Теперь мы имеем пару нетривиальных делителей числа n - q и n/q
		//Разложим каждый из них на множители (углубим рекурсию)
		//множитель q
		pollard_factor(q,tmp_rez,fact_num);
		for (i=0;i<fact_num;i++)
			results[rez_count++]=tmp_rez[i];
		//множитель n/q
		pollard_factor(n/q,tmp_rez,fact_num);
		for (i=0;i<fact_num;i++)
			results[rez_count++]=tmp_rez[i];
		}
}

//Считая, что факторизованное число было порядком группы
//По разложению на простые числа найдём возможные порядки подгрупп
void find_subgroup_orders(ZZ * input_vec, long input_len, ZZ * output_vec, long & output_len) {
ZZ * prime_vec = new ZZ [input_len];
long prime_ub;
long i,j;
bool flag;
for (i=0;i<input_len;i++) {
	flag = false;
	for (j=0;j<prime_ub;j++) {
		if (prime_vec[j]==input_vec[i]) {
			flag = true;
			output_vec[j] *= input_vec[i];
			break;
			}
		}
	if (!flag) {
		prime_vec[prime_ub]=input_vec[i];
		output_vec[prime_ub]=input_vec[i];
		prime_ub++;
		}
	}
output_len = prime_ub;
delete [] prime_vec;
}


/*
int main() {
ifstream fin("input.txt");
ofstream fout("output.txt");
ZZ n;
ZZ x0,q;

fin >> n;

ZZ results[MAX_FACT];
long rez_count,i;

pollard_factor(n,results,rez_count);
ZZ subgroups[MAX_FACT];
long subgroups_count;
find_subgroup_orders(results,rez_count,subgroups,subgroups_count);

for (i=0;i<subgroups_count;i++) {
	fout << subgroups[i] << endl;
	}
fin.close();
fout.close();
}
*/
