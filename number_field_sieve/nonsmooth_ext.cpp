#include "nonsmooth_ext.h"
#include "sieve.h"
#include <fstream>
using namespace std;

NTL_CLIENT
//Неплохо бы переписать под параллельные вычисления
bool make_smooth(ZZ x, ZZ p, ZZ smooth_checker, ZZ * basis, long basis_size, //in
		 int crank, int csize, int block_size,  //in
		ZZ& rez, ZZ& rez_power, ZZ * result_powers) //out
{
	rez_power=crank+1;
	rez = PowerMod(x,rez_power,p);
	ZZ mul = PowerMod(x,csize,p);
	int found,found_sum;
	int i;

	do {
	  found=0;
	  for (i=0;i<block_size;i++) {
		if (sieving(rez,p,smooth_checker,basis,basis_size,result_powers)) {
			found=1;
			break;
			}
		rez = MulMod(rez,mul,p);
		rez_power+=csize;
		}
	MPI_Allreduce(&found,&found_sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	} while (found_sum==0);

	int t,who_found=0;
	MPI_Status status;
      if (crank==0) {
	for (i=1;i<csize;i++) {
	  MPI_Recv(&t,1,MPI_INT,i,1000,MPI_COMM_WORLD,&status);
	  if (t!=0)
	    who_found=i;
	  }
	}
      else {
	MPI_Send(&found,1,MPI_INT,0,1000,MPI_COMM_WORLD);
	}

MPI_Bcast(&who_found,1,MPI_INT,0,MPI_COMM_WORLD);

if (crank==0) {
  if (who_found==0)
      return true;
  else {
      MPI_NTL_Recv(rez,who_found);
      MPI_NTL_Recv(rez_power,who_found);
      for (i=0;i<basis_size;i++)
	MPI_NTL_Recv(result_powers[i],who_found);
      } 
  }
else {
  if (crank==who_found) {
      MPI_NTL_Send(rez,0);
      MPI_NTL_Send(rez_power,0);
      for (i=0;i<basis_size;i++)
	MPI_NTL_Send(result_powers[i],0);
    }
    
  }

}

bool solve_matrix(ZZ ** base, int lines, int columns, ZZ r, //in
		  ZZ& alpha, ZZ& beta) { //out
int i,j,k;

//Получили матрицу. Первые две строки - коэффициенты разложения a,b
//Остальные - вектора разложения 1
//Проверим на наличие вырожденных столбцов
bool degraded;
for (j=0;j<columns;j++) {
  degraded=true;
  for (i=2;i<lines;i++)
    if (base[i][j]!=0) 
      degraded=false;
  if (degraded) {
    cout << "В матрице new_base найден выродившийся столбец" << endl;
    return false;
    }
  }


//Транспонируем матрицу
ZZ ** t_base = new ZZ*[columns];

for (j=0;j<columns;j++)
  t_base[j] = new ZZ[lines];

for (i=0;i<lines;i++)
  for (j=0;j<columns;j++) {
    t_base[j][i]=base[i][j];
    }
k=lines;
lines=columns;
columns=k;
//if (columns>(lines+1)) 
//  columns=lines+1;
//Завершили транспонирование

bool * processed = new bool[lines];
for (i=0;i<lines;i++)
  processed[i]=false;

int base_line;
ZZ tz;
for (j=2;j<columns;j++) { //Устраняем все лишние столбцы (кроме первых двух)
  //Выбираем опорную строку
  base_line=-1;
  for (i=0;i<lines;i++) {
    if (!processed[i] && t_base[i][j]!=0) {
      base_line=i;
     break;
      }
    }

  if (base_line==-1) continue;

  processed[base_line]=true;

  for (i=0;i<lines;i++) {
    if (!processed[i]) {
      tz=t_base[i][j];
      for (k=0;k<columns;k++) {
	t_base[i][k]=tz*t_base[base_line][k] - t_base[base_line][j]*t_base[i][k];
	t_base[i][k]%=r;
	}
      }
    }
  }

alpha=0;
beta=0;
for (i=0;i<lines;i++) {
   if (!processed[i]) {
    alpha+=t_base[i][0];
    beta+=t_base[i][1];
    }
   if (alpha!=0 && beta!=0) break;
   }
for (j=0;j<lines;j++)
  delete[] t_base[j];
delete[] t_base;
delete[] processed;

if (alpha==0 || beta==0) return false;

beta=NegateMod(beta,r);
tz=alpha;
alpha=beta;
beta=tz;

return true;
}


#ifdef MODULE_TEST
int main() {

long d1_size=5000;

ifstream ftask("task.txt");
ZZ a,b,p,r;
ftask >> a >> b >> p >> r;
ftask.close();

//Построение базиса D1
ZZ * d1 = new ZZ[d1_size];
ZZ d1_smooth_checker=to_ZZ(1);

d1[0]=to_ZZ(-1); //Включая число -1, ибо так сказал Великий Ростовцев
for (long i=1;i<d1_size;i++) {
	NextPrime(d1[i], d1[i-1]+1);
	d1_smooth_checker *= power(d1[i],(long)ceil(log(p)/log(d1[i])));
	}

ZZ rez_a,rez_power_a, *result_powers_a=new ZZ[d1_size];
ZZ rez_b,rez_power_b, *result_powers_b=new ZZ[d1_size];

make_smooth(a,p,d1_smooth_checker, d1, d1_size,  //in
		rez_a, rez_power_a, result_powers_a); //out

cout << rez_a << " " << rez_power_a << endl;

make_smooth(b,p,d1_smooth_checker, d1, d1_size,  //in
		rez_b, rez_power_b, result_powers_b); //out

cout << rez_b << " " << rez_power_b << endl;

return 0;
}
#endif
