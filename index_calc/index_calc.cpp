/*
Распределёная программа, осуществляющая поиск дискретного логарифма в GF(p)
методом базы разложения.
Автор: Игорь Сидоров, аспирант каф. БИТ ТТИ ЮФУ
*/

#include "mpi.h" 

#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include "sieve.h"
#include "gauss.h"
#include "mpi_ntl_routines.h"

NTL_CLIENT

#define MES_NN_PRIMES 1
#define MES_WHICH_COL 2
#define MES_WHO_SEND 3
#define MES_FOUND_B_DECOMP 4

int main(int argc, char ** argv) {
//основание, степень, модуль, порядок группы
ZZ a,b,p,r; 
//Базис разложения
ZZ * basis;
//Размер базиса
long basis_size;
//для ускорения решета
ZZ basis_mul;

//Индексы циклов
long i,j,k,i2,j2;
//Отсечки времени
double t1,t2;

//Номер процесса, число процессов
int crank,csize;

//Число дополнительных строк матрицы.
//Расширенная матрица имеет размер (basis_size+1 столбцов) * (basis_size+matrix_addon_size строк)
int matrix_addon_size;
int start_power;

int a_max_count,b_max_count;
int save_matrix;

bool * nn_primes;
bool * root_nn_primes;

MPI_Status status;

MPI_Init( &argc, &argv );


ifstream ftask("task.txt");
ifstream fparam("param.txt");

//Вводим исходные данные
ftask >> a >> b >> p >> r;

//Считываем параметры
fparam >> basis_size;
fparam >> matrix_addon_size;
fparam >> start_power;
fparam >> a_max_count >> b_max_count;
fparam >> save_matrix;
ftask.close();
fparam.close();

//Строим базис
basis = new ZZ[basis_size];

basis_mul=to_ZZ(1);
basis[0]=to_ZZ(-1); //Включая число -1, ибо так сказал Великий Ростовцев
for (i=1;i<basis_size;i++) {
	NextPrime(basis[i], basis[i-1]+1);
	if (i==1) basis_mul *= power(basis[i],30);
	else if (i==2) basis_mul *= power(basis[i],20);
	else if (i==3) basis_mul *= power(basis[i],15);	
	else if (i==4) basis_mul *= power(basis[i],10);	
	else if (i==5) basis_mul *= power(basis[i],10);	
	else if (i==6) basis_mul *= power(basis[i],10);	
	else if (i==7) basis_mul *= power(basis[i],10);		
	else if (i==8) basis_mul *= power(basis[i],8);		
	else if (i==9) basis_mul *= power(basis[i],8);		
	else if (i==10) basis_mul *= power(basis[i],8);		
	else if (i==11) basis_mul *= power(basis[i],5);	
	else           basis_mul *= power(basis[i],3);	
	}

//basis[0]=to_ZZ(-1);
//for (i=1;i<basis_size;i++) {
//	basis_mul *= power(basis[i],(long)ceil(log(p)/log(basis[i])));
//	}

//Определяем какое наше место
MPI_Comm_rank( MPI_COMM_WORLD, &crank );
//и сколько нас
MPI_Comm_size( MPI_COMM_WORLD, &csize );

//сколько строк матрицы будет в нашем стрипе
int matrix_size=basis_size + matrix_addon_size;

//База разложения (стрип большой злой матрицы)
ZZ ** base = new ZZ * [matrix_size];

//Текущий размер нашего стрипа
int base_ub=0;

//Строка этой матрицы (для функции просева через решето)
ZZ * result_powers = new ZZ [basis_size+1];

/* Текущий кандидат для просеивания через решето, множитель, показатель текущего кандидата
cur_pow = a^x mod p
*/
ZZ cur_pow,multiplier,x;

//Инициализация текущего кандидата
x=to_ZZ(crank+start_power);
cur_pow=PowerMod(a,x,p);

multiplier=PowerMod(a,csize,p);

//Цикл просеивания через решето
//Пока не заполнен стрип
if (crank==0)
	cout << "NumBytes(Basis_mul)= " << NumBytes(basis_mul) << endl;

t1=MPI_Wtime();

int a_count;
int base_ub_sum=0;
while (base_ub_sum<matrix_size) {
	for (a_count=0;(a_count<a_max_count) && (base_ub<matrix_size);a_count++) {
		//Проверяем текущего кандидата на решете
		if (sieving(cur_pow, p, basis_mul, basis,basis_size,result_powers)) {
			//Если он B-гладкий
			//Добавим в наш стрип матрицы
			base[base_ub]= new ZZ [basis_size+1];
			for (i=0;i<basis_size;i++)
				base[base_ub][i]=result_powers[i];
			//  также в матрицу пойдёт показатель кандидата
			// (как столбец свободных членов)
			base[base_ub][basis_size]=x;
			base_ub++;
			};
		//Переходим к следующему кандидату	
		cur_pow=MulMod(cur_pow,multiplier,p);
		x+=to_ZZ(csize);
		}
	//Считаем, сколько гладких степеней нашли все процессы
	MPI_Allreduce(&base_ub,&base_ub_sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	}

t2=MPI_Wtime();
cout << "Процесс #" << crank << " из " << csize << " просеял участок матрицы размером " << base_ub << " за " << t2-t1 << " секунд" << endl;

//Ожидаем, пока другие процессы не заполнят свои стрипы
MPI_Barrier(MPI_COMM_WORLD);


nn_primes = new bool [basis_size];
root_nn_primes = new bool [basis_size];

for (j=0;j<basis_size;j++) {
	nn_primes[j] = false;
	root_nn_primes[j] = false;
	}


for (i=0;i<base_ub;i++)
	for (j=0;j<basis_size;j++)
		if (base[i][j]!=0)
			nn_primes[j]=true;

int nn_size=0;

//Согласовываем переменные, которые будут участвовать в рассмотрении
if (crank==0) {
	for (j=0;j<basis_size;j++)
		root_nn_primes[j]=nn_primes[j];
	for (i=1;i<csize;i++) {
		MPI_Recv(nn_primes,basis_size*sizeof(bool),MPI_CHAR,i,MES_NN_PRIMES,MPI_COMM_WORLD,&status);
		for (j=0;j<basis_size;j++)
			root_nn_primes[j]|=nn_primes[j];
		}
	for (j=0;j<basis_size;j++)
		if (root_nn_primes[j]) nn_size++;
	cout << "Получил вектор ROOT_NN_PRIMES. Ненулевых компонентов: " << nn_size << endl;
	}
else {
	MPI_Send(nn_primes,basis_size*sizeof(bool),MPI_CHAR,0,MES_NN_PRIMES,MPI_COMM_WORLD);
	}

MPI_Bcast(root_nn_primes,basis_size*sizeof(bool),MPI_CHAR,0,MPI_COMM_WORLD);

//Просеиваем степени числа B
//Распределённый вариант
if (crank==0) {
	t1=MPI_Wtime();
	}

ZZ b_power=to_ZZ(crank);
ZZ b_cur_pow=PowerMod(b,b_power,p);
ZZ b_mul = PowerMod(b,csize,p);
ZZ * b_base = new ZZ [basis_size+1];

bool b_flag=false;
int who_found_b=1,am_i;
int b_count;

who_found_b=-1;

while (who_found_b==-1) {
	b_count=0;

	//Собственно просеивание числа b
	while ((!b_flag) && (b_count<b_max_count)) {
		b_count++;
		MulMod(b_cur_pow,b_cur_pow,b_mul,p);
		b_power+=csize;
		if (sieving(b_cur_pow, p, basis_mul, basis,basis_size,result_powers)) {
			b_flag=true;
			for (j=0;j<basis_size;j++)
				if (result_powers[j]!=0 && root_nn_primes[j]==false)
					b_flag=false;
			}
		}
	
	//Нашли число - готовимся поделиться с остальными
	if (b_flag) {
		for (j=0;j<basis_size;j++)
			b_base[j]=result_powers[j];
		b_base[basis_size]=0;
		}
	
	who_found_b=-1;
	if (crank==0) {
		if (b_flag) 
			who_found_b=0;
		for (i2=1;i2<csize;i2++) {
			MPI_Recv(&am_i,1,MPI_INT,i2,MES_FOUND_B_DECOMP,MPI_COMM_WORLD,&status);
			if (am_i!=-1 && who_found_b==-1)
				who_found_b=i2;
			}
		}
	else {
		am_i=b_flag?crank:-1;
		MPI_Send(&am_i,1,MPI_INT,0,MES_FOUND_B_DECOMP,MPI_COMM_WORLD);
		}
	//Рассылаем что нам известно о том, нашёл ли кто разложение b
	MPI_Bcast(&who_found_b,1,MPI_INT,0,MPI_COMM_WORLD);
	}

//Выдаём всем процессам найденное разложение B

/* Старый вариант
for (i2=0;i2<=basis_size;i2++)
	MPI_NTL_Bcast(b_base[i2],who_found_b);

MPI_NTL_Bcast(b_power,who_found_b);
*/

unsigned char * temp_buf;
int temp_buf_size,temp_buf_ub;

temp_buf = new unsigned char [(8+NumBytes(p)+2)*(basis_size+2)];
temp_buf_ub=0;
for (i2=0;i2<=basis_size;i2++)
  NTL_Pack(b_base[i2],temp_buf,temp_buf_ub);
NTL_Pack(b_power,temp_buf,temp_buf_ub);
temp_buf_size=temp_buf_ub;
//Широковещательно рассылаем
Buf_bcast(temp_buf,temp_buf_size,who_found_b);
//Распаковываем
temp_buf_ub=0;
for (i2=0;i2<=basis_size;i2++)
	NTL_Unpack(b_base[i2],temp_buf,temp_buf_ub,temp_buf_size);
NTL_Unpack(b_power,temp_buf,temp_buf_ub,temp_buf_size);
delete [] temp_buf;

if (crank==0) {
	t2=MPI_Wtime();
	cout << "Процесс # " << who_found_b << " нашёл разложение b за время " << t2-t1 << " секунд" << endl;
	}


ofstream fmatrix_output;
char file_name[512];
if (save_matrix) {
	sprintf(file_name,"strip_%ld",crank);
	fmatrix_output.open(file_name);
		for (i=0;i<base_ub;i++) {
			for (j=0;j<basis_size+1;j++)
				fmatrix_output << base[i][j] << " ";
			fmatrix_output << endl;
			}
	fmatrix_output.close();
	if (crank==0) {
		fmatrix_output.open("b_decomp");
		for (i=0;i<basis_size;i++)
			fmatrix_output << b_base[i] << " ";
		fmatrix_output << b_power << endl;
		fmatrix_output.close();
		}
	}


//Вторая часть решения - цикл Гауссова исключения
if (crank==0) 
	t1 = MPI_Wtime();
//Текущий результат поиска строки, в которой заданный элемент ненулевой
int cur_comp;
//Сохранённый про запас результат поисков
int back_cur_comp;
//Какой процесс рассылает опорную строку. Если -1 - рассылать некому - пропускаем данную переменную
int who_send;

bool * processed_vectors = new bool [base_ub];
for (i=0;i<base_ub;i++)
	processed_vectors[i]=false;

//Опорная строка, по которой процессы будут обрабатывать свои стрипы
ZZ * active_string = new ZZ[basis_size+1];

for (j=0;j<basis_size;j++) { //Цикл по переменным

	if (root_nn_primes[j]==false) continue;

	//Находим в каждом стрипе вектор с ненулевым заданным (j) компонентом
	//Если не нашли, вернём -1
	cur_comp=-1;
	who_send=-1;
	for (i=0;i<base_ub;i++) {
		if (processed_vectors[i]==false) {

			if (base[i][j]!=0) {
				cur_comp = i;
				break;
				}
			}
		}

	back_cur_comp = cur_comp;

	if (crank!=0) {
		//Если мы не корневой процесс - пошлём результаты наших поисков
		MPI_Send(&cur_comp,1,MPI_INT,0,MES_WHICH_COL,MPI_COMM_WORLD);
		}
	else {
		//А если корневой - объединим усилия всех процессов
		k=0;
		if ((cur_comp!=-1) && (who_send==-1)) who_send=k;

		for (k=1;k<csize;k++) {
			MPI_Recv(&cur_comp,1,MPI_INT,k,MES_WHICH_COL,MPI_COMM_WORLD,&status);
			if (cur_comp!=-1 && who_send==-1) who_send=k;
			}
		}

	//Разошлём всем процессам информацию о том, кто будет рассылать вектор
	MPI_Bcast(&who_send,1,MPI_INT,0,MPI_COMM_WORLD);

	//Если некому рассылать
	if (who_send==-1) {
		continue;
		}

	//Если обрабатываемый вектор принадлежит текущему процессу
	if (who_send==crank) {
		//Помечаем вектор как обработанный
		processed_vectors[back_cur_comp]=true;
		//Копируем его
		for (k=0;k<=basis_size;k++)
			active_string[k]=base[back_cur_comp][k];
		}


	//Рассылаем всем активный вектор
	for (k=0;k<j;k++)
		active_string[k]=0;

/* Старый вариант
	for (k=j;k<=basis_size;k++)
		MPI_NTL_Bcast(active_string[k],who_send);
*/

temp_buf = new unsigned char [(8+NumBytes(p)+2)*(basis_size+1)];
temp_buf_ub=0;
for (k=j;k<=basis_size;k++)
  NTL_Pack(active_string[k],temp_buf,temp_buf_ub);
temp_buf_size=temp_buf_ub;
//Широковещательно рассылаем
Buf_bcast(temp_buf,temp_buf_size,who_send);
//Распаковываем
temp_buf_ub=0;
for (k=j;k<=basis_size;k++)
	NTL_Unpack(active_string[k],temp_buf,temp_buf_ub,temp_buf_size);
delete [] temp_buf;

	//Теперь для каждого необработанного вектора стрипа матрицы
	//Выполняем Гауссово преобразование
	for (i2=0;i2<base_ub;i2++) {
		if (processed_vectors[i2]==false) {
			make_vector(base[i2],active_string,j,basis_size+1,r);
			}
		}

	//Такому же преобразованию подвергается вектор B
	if (crank == 0) {
		if (b_base[j]!=0) {
			make_vector(b_base,active_string,j,basis_size+1,r);
			b_power=MulMod(b_power,active_string[j],r);
			}

		}
	}

//Выводим матрицу !после! Гауссова исключения
//Чтобы посчитать кол-во нулевых (лин. зависимых) строк
	bool zero_flag;
	int zero_count=0,sum_zero_count;
	for (i=0;i<base_ub;i++) {
		zero_flag=true;
		for (j=0;j<basis_size+1;j++) {
			if (base[i][j]!=0) zero_flag=false;
			}
		if (zero_flag) zero_count++;
	    }
MPI_Reduce(&zero_count,&sum_zero_count,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
if (crank==0) {
    cout << "Число выродившихся строк: " << sum_zero_count << endl;
    }

ZZ ans,divisor,a_power,corr,t_norm;
int test;
if (crank == 0 ) {
	t2=MPI_Wtime();
	cout << "Работа с матрицей заняла " << t2-t1 << " секунд" << endl;

	//Третий этап - решение уравнения
	//Проверим, все ли компоненты разложения степени b по базису равны нулю. 
	test=0;
	for (i2=0;i2<basis_size;i2++)
		if (b_base[i2]!=0) {
			//cout << "B component #" << i2 << " is nonzero!" << endl;
			test++;
			}

	//Если нет - при решении произшла ошибка
	if (test!=0) {
		cout << "Найдено " << test <<" ненулевых компонент в векторе B!" << endl;
		cout << "Решение невозможно найти!" << endl;
		MPI_Finalize();
		exit(1);
		}

	a_power = b_base[basis_size];
	cout << "a_power = " << a_power << " b_power = " << b_power << endl;
	//Ответ находится из соотношения a^a_power * b^b_power = 1 (mod p)
	//Отсюда a ^ (a_power + x*b_power) = 1 (mod p)
	//(a_power + x*b_power) = 0 (mod r)
	// x = b_power^-1 - a_power (mod r)
	//Эти проверки позволяют корректно решать задачу, если r - составное
	//Обычно a - образующая, и r=(p-1)=2*p', где p' - простое
	divisor = GCD(b_power,r);
	cout << "GCD(b_power,r) = " << divisor << endl;
	r/=divisor;
	//Правило гласит, что a=b (mod rs), тогда и только тогда, когда a=b (mod r) и a=b (mod s)
	//При этом НОД(r,s)=1, иначе правило работает не всегда
	do {
		t_norm=GCD(divisor,r);
		r/=t_norm;
		divisor*=t_norm;
		} while (t_norm!=1);

	cout << "Делитель r = " << divisor << endl;
	ans = MulMod(InvMod(b_power%r,r),NegateMod(a_power%r,r),r);

	//А теперь скорректируем ответ, если divisor!=1
	for (corr=0;corr<(r*divisor);corr+=r) {
		if (PowerMod(a,ans+corr,p)==b) {
			ans+=corr;
			break;
			}
		}
	cout << "Ответ = " << ans << endl;

	//Проверим, правильный ли мы получили ответ
	if (PowerMod(a,ans,p)==b) {
		cout << "Ответ верен" << endl;
		}
	else {
		cout << "Ответ НЕВЕРЕН!" << endl;
		}

	}
//Корректное завершение работы программы
MPI_Finalize();
return 0;
}
