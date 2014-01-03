/*
Распределёная программа, осуществляющая поиск дискретного логарифма в GF(p)
методом решета числового поля
Автор: Игорь Сидоров, аспирант каф. БИТ ТТИ ЮФУ
*/

#include "mpi.h"
#include "sieve.h"
#include "gauss.h"
#include "mpi_ntl_routines.h"
#include "nonsmooth_ext.h"
#include <NTL/ZZ.h>
#include <iostream>
#include <fstream>
using namespace std;

NTL_CLIENT 

#define MES_WHICH_COL 2

#include "algebraic_number.h"

/** Функция генерации алгебраического числа по заданному номеру
@param id - порядковый номер вектора показателей
@param d2 - базис алгебраических чисел D2
@param d2_size - размер базиса
@param bound - максимальная степень
@param result_powers - массив для возвращения результата (вектора показателей)
@return алгебраическое число (элементы базиса D2, взятые с соответствующими степенями)
*/
algebraic_number convert_id_to_alg_number(ZZ id, algebraic_number * d2, int d2_size, ZZ bound, ZZ * result_powers) {
	int i;
	result_powers[0] = RandomBnd(2);
	result_powers[1] = RandomBnd(2);

	for (i=2;i<d2_size;i+=2) {
		if (RandomBnd(2)==0) {
			//Выбираем число или сопряжённое к нему, но не оба сразу
			result_powers[i]=RandomBnd(bound);
			result_powers[i+1]=0;
			}
		else {
			result_powers[i]=0;
			result_powers[i+1]=RandomBnd(bound);
			}
		}
	//Получение алгебраического числа по показателям элементов базиса d2
	algebraic_number result(1,0);
	for (i=0;i<d2_size;i++)
		result = result * d2[i].power(result_powers[i]);
	return result;
	}

/** Зануление указанного индекса с помощью метода Гаусса
@param ind Индекс (столбец), который мы зануляем
@param base Матрица показателей
@param base_ub Число строк матрицы
@param var_count Число столбцов матрицы
@param processed_vectors Булев вектор, в котором помечается, участвовала ли строка в качестве опорной
@param crank Ранг (номер) процесса
@param csize Число процессов
@param r Порядок элемента. 0 - вычисления без модуля.
*/
void make_zero_index(long ind, ZZ ** base, long base_ub, long var_count, bool * processed_vectors, long crank, long csize, ZZ r) {
	int i,k;
	MPI_Status status;

	//Текущий результат поиска строки, в которой заданный элемент ненулевой
	int cur_comp;

	//Сохранённый про запас результат поисков
	int back_cur_comp;

	//Какой процесс рассылает опорную строку.
	//Если -1 - рассылать некому - пропускаем данную переменную
	//Это может привести к проблемам в методе базы разложения
	//но не критично в решете числового поля
	int who_send;

	//Опорная строка, по которой процессы будут обрабатывать свои стрипы
	ZZ * active_string = new ZZ[var_count];

	//Находим в каждом стрипе вектор с ненулевым заданным (ind) компонентом
	//Если не нашли, вернём -1
	cur_comp=-1;
	who_send=-1;
	for (i=0;i<base_ub;i++) {
		if (processed_vectors[i]==false) {
			if (base[i][ind]!=0) {
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
//		if (crank==0)
//			cout << "Компонент #" << ind << " невозможно занулить!" << endl;
		if (active_string!=NULL)
			delete [] active_string;
		return;
		}

	//Если обрабатываемый вектор принадлежит текущему процессу
	if (who_send==crank) {
		//Помечаем вектор как обработанный
		processed_vectors[back_cur_comp]=true;
		//Копируем его
		for (k=0;k<var_count;k++)
			active_string[k]=base[back_cur_comp][k];
		}


	//Рассылаем всем активный вектор
	/* Старый вариант - рассылаем одиночными сообщениями
	for (k=0;k<var_count;k++)
		MPI_NTL_Bcast(active_string[k],who_send);
	*/

	//Новый вариант - с буферизацией
	int buf_reserved=(8+NumBytes(r)+2)*var_count;
	unsigned char *temp_buf=new unsigned char[buf_reserved];
	int temp_buf_size,temp_buf_ub=0;
	//Упаковываем
	for (k=0;k<var_count;k++)
		NTL_Pack(active_string[k],temp_buf,temp_buf_ub);

	temp_buf_size=temp_buf_ub;
	//Широковещательно рассылаем
	Buf_bcast(temp_buf,temp_buf_size,who_send);
	//Распаковываем
	temp_buf_ub=0;
	for (k=0;k<var_count;k++)
		NTL_Unpack(active_string[k],temp_buf,temp_buf_ub,temp_buf_size);

        delete[] temp_buf;
	//Конец кода рассылки опорной строки

	//Теперь для каждого необработанного вектора стрипа матрицы
	//Выполняем Гауссово преобразование
	int i2;
	for (i2=0;i2<base_ub;i2++) {
		if (processed_vectors[i2]==false) {
			make_common_vector(base[i2],active_string,ind,var_count,r);
			}
		}

	if (active_string!=NULL)
		delete [] active_string;
};

//Главная функция программы
int main(int argc, char ** argv) {

//Инициализируем MPI
MPI_Init(&argc, &argv);
//Номер процесса, число процессов
int crank,csize;
long i,j,k;

MPI_Comm_rank(MPI_COMM_WORLD,&crank);
MPI_Comm_size(MPI_COMM_WORLD,&csize);

SetSeed(to_ZZ(crank));

//Шаг 0 - подготовительные операции
//Ввод исходных данных
//основание, степень, модуль, порядок группы
ZZ a,b,p,r; 
ifstream ftask("task.txt");
//Вводим исходные данные
ftask >> a >> b >> p >> r;
ftask.close();

//Ввод параметров алгоритма
ZZ sqr_alpha,m;
ZZ d2_index_bound;
long d1_size,d2_size;
long matrix_addon_count;
long cand_max_count;

ifstream fparam("param.txt");
fparam >> d1_size >> d2_size;
fparam >> matrix_addon_count;
fparam >> d2_index_bound;
fparam >> cand_max_count;
fparam.close();

//Построение базиса D1
ZZ * d1 = new ZZ[d1_size];
ZZ d1_smooth_checker=to_ZZ(1);

d1[0]=to_ZZ(-1); //Включая число -1, ибо так сказал Великий Ростовцев
for (i=1;i<d1_size;i++) {
	NextPrime(d1[i], d1[i-1]+1);
	if (i==1) d1_smooth_checker *= power(d1[i],30);
	else if (i==2) d1_smooth_checker *= power(d1[i],20);
	else if (i==3) d1_smooth_checker *= power(d1[i],15);	
	else if (i==4) d1_smooth_checker *= power(d1[i],10);	
	else if (i==5) d1_smooth_checker *= power(d1[i],10);	
	else if (i==6) d1_smooth_checker *= power(d1[i],10);	
	else if (i==7) d1_smooth_checker *= power(d1[i],10);		
	else if (i==8) d1_smooth_checker *= power(d1[i],8);		
	else if (i==9) d1_smooth_checker *= power(d1[i],8);		
	else if (i==10) d1_smooth_checker *= power(d1[i],8);		
	else if (i==11) d1_smooth_checker *= power(d1[i],5);	
	else           d1_smooth_checker *= power(d1[i],3);	
	}

//d1[0]=to_ZZ(-1);
//for (i=1;i<d1_size;i++) {
//	d1_smooth_checker *= power(d1[i],(long)ceil(log(p)/log(d1[i])));
//	}

if (crank==0)
    cout << "NumBytes(d1_smooth_checker)=" << NumBytes(d1_smooth_checker) << endl;

sqr_alpha=0;

for (i=1;i<d1_size;i++) {
	if (Jacobi(-d1[i]%p,p)==1) {
		sqr_alpha=-d1[i];
		break;
		}
	}

if (sqr_alpha==0) {
	cout << "Не могу построить расширенное поле или установить изоморфизм" << endl;
	MPI_Finalize();
	exit(0);
	}

//Инициализируем класс алгебраических чисел
algebraic_number::init_sqr_alpha(sqr_alpha);

//Найдём m
m = SqrRootMod(sqr_alpha%p,p);
if (crank==0) {
	cout << "sqr_alpha = " << sqr_alpha << endl;
	cout << "m = " << m << endl;
	}

//Построение базиса D2
algebraic_number * d2 = new algebraic_number[d2_size+1];
long d2_ub=0;
bool d2_exit=false;
ZZ tc,td;
algebraic_number d2_cand;

//Включая идеал -1, ибо так подразумевал Великий Ростовцев
d2[d2_ub++]=to_algebraic_number(to_ZZ(-1),to_ZZ(0));

d2[d2_ub++]=to_algebraic_number(to_ZZ(0),to_ZZ(1));

for (i=1;d2_ub<d2_size;i++) {
	for (j=1;(j<i) && (d2_ub<d2_size);j++) {
		tc=j;
		td=i-j;
		if (td==0) continue;
		d2_cand = to_algebraic_number(tc,td);
		if (ProbPrime(d2_cand.norm())) { //Идеал будет простым, если его норма - простое число
			d2[d2_ub]=d2_cand;
			d2_ub++;
			d2[d2_ub++]=d2_cand.get_mate();
			}
		}
	}

d2_size = d2_ub;

//Сделаем a,b D1-гладкими
ZZ a_smooth,a_smooth_ind, *a_smooth_decomp = new ZZ[d1_size];
ZZ b_smooth,b_smooth_ind, *b_smooth_decomp=new ZZ[d1_size];

bool *d1_must_leave=new bool[d1_size];
int d1_must_leave_nonzero=0;
double t1,t2; //Отсечки времени
if (crank==0) 
  t1=MPI_Wtime();

  //Сделаем гладким число a
  make_smooth(a, p, d1_smooth_checker, d1, d1_size, crank, csize, cand_max_count, //in
		a_smooth, a_smooth_ind, a_smooth_decomp); //out
  //Сделаем гладким число b
  make_smooth(b, p, d1_smooth_checker, d1, d1_size, crank, csize, cand_max_count, //in
		b_smooth, b_smooth_ind, b_smooth_decomp); //out

  if (crank==0) {
  //Сформируем вектор степеней, которые мы не будем исключать
  for (i=0;i<d1_size;i++) {
    d1_must_leave[i]=(a_smooth_decomp[i]!=0) || (b_smooth_decomp[i]!=0);
    if (d1_must_leave[i]) {
      cout << d1[i] << " ";
      d1_must_leave_nonzero++;
      }
    }
  t2=MPI_Wtime();
  cout << endl;
  cout << "Шаг построения гладкой задачи занял " << t2-t1 << " секунд" << endl;
  cout << "Число ненулевых показателей в разложениях: " << d1_must_leave_nonzero << endl;
  }

MPI_Bcast(d1_must_leave,d1_size,MPI_CHAR,0,MPI_COMM_WORLD);
MPI_Bcast(&d1_must_leave_nonzero,1,MPI_INT,0,MPI_COMM_WORLD);

//Шаг 1 - построение пар (c,d), таких, что
//c+dm (mod p) - D1 - гладкое
//с+d*alpha - D2 - гладкое

ZZ c,d;
//Базы разложения
//Первый индекс - номер разложенного эл-та, второй - номер простого числа (идеала) в базисах D1, D2

//Стрип матрицы
ZZ **base;

//Общий размер матрицы
int base_size = d1_size + d2_size + 2 + matrix_addon_count; //#D1 + #D2 + n

base = new ZZ * [base_size]; //Выделяем место, пока только под указатели

int base_ub=0;  //Текущее число элементов в стрипе данного процесса

int cand_count;  //Индекс по числу кандидатов между проверками размера матрицы

int base_ub_sum=0; //Число элементов во всех стрипах (по всем процессам)

ZZ x=to_ZZ(crank+2);

//Разложение элементов-кандидатов по базисам D1,D2
ZZ * d1_powers = new ZZ[d1_size];
ZZ * d2_powers = new ZZ[d2_size];

ZZ z_cand;
algebraic_number alg_cand;

t1=MPI_Wtime();

while (base_ub_sum<base_size) { //Пока не заполним всю матрицу
	//Проверяем заданное число кандидатов
	for (cand_count=0;(cand_count<cand_max_count) && (base_ub<base_size);cand_count++) { 
		//Пара (c,d) в виде алгебраического числа c+d*alpha 
		//Сгенерирована d2-гладкой
		alg_cand = convert_id_to_alg_number(x,d2,d2_size,d2_index_bound,d2_powers);
		if (((alg_cand.get_c()==0) || (alg_cand.get_c()==1)) && (alg_cand.get_d()==0)) {
			x+=to_ZZ(csize);
			continue;
			}

		//Получим целое число из GF(p)
		z_cand = (alg_cand.get_c() + alg_cand.get_d()*m)%p;
		//Проверим его на d1-гладкость
		if (sieving(z_cand, p, d1_smooth_checker, d1,d1_size,d1_powers)) {
			//Если в разложении элементов присутствует -1
			if ((d1_powers[0]==1) && (d2_powers[0]==1)) {
				x+=to_ZZ(csize);
				continue;
				}
			//Если гладкое - выделим память под новую строку матрицы
			base[base_ub]=new ZZ[d1_size+d2_size];
			//Скопируем показатели элементов базиса d1,d2
			for (i=0;i<d1_size;i++)
				base[base_ub][i] = d1_powers[i];
			for (i=0;i<d2_size;i++)
				base[base_ub][i+d1_size] = d2_powers[i];
			//Увеличим число сохранённых элементов
			base_ub++;
			};
		//Переходим к следующему кандидату
		x+=to_ZZ(csize);
		}
	//Считаем, сколько гладких элементов нашли все процессы
	MPI_Allreduce(&base_ub,&base_ub_sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	}

t2=MPI_Wtime();
cout << "Процесс #" << crank << " из " << csize << " нашел полосу размером " << base_ub << " строк за  " << t2-t1 << " секунд" << endl;

//Шаг 2 - линейная алгебра
bool * processed_vectors = new bool [base_ub];

for (i=0;i<base_ub;i++)
	processed_vectors[i]=false;

//Занулим показатели при элементах базиса d2
//Наша задача получить степень алгебраического числа, кратную r
//А если мы загоним все элементы базиса d2 в 0, то 1=1^0, 0 кратен r
if (crank == 0)
    cout << "Начал исключение элементов базиса D2" << endl;
//цикл Гауссова исключения
long ind;
if (crank==0) 
	t1 = MPI_Wtime();
for (ind=d1_size;ind<d1_size+d2_size;ind++) { //Цикл по переменным
	make_zero_index(ind, base, base_ub, d1_size+d2_size, processed_vectors, crank, csize,r);
	}

//Занулим показатели при элементах базиса d1
//кроме тех, которые входят в разложение a и b

//Собственно зануление
if (crank == 0)
    cout << "Начал исключение элементов базиса D1" << endl;
for (ind=0;ind<d1_size;ind++) {
	if (!d1_must_leave[ind]) {
		make_zero_index(ind, base, base_ub, d1_size+d2_size, processed_vectors, crank, csize,r);
		}
	}

if (crank==0) {
  t2=MPI_Wtime();
  cout << "Обработка матрицы заняла " << t2-t1 << " секунд." << endl;
  }


//Получим коэффициенты для составления сравнения
//s-это коэффициент при a, t-коэффициент при b
//Выражение имеет вид (a^s)*(b^t)=1(mod p)

ZZ my_s=to_ZZ(0), my_t=to_ZZ(0);
ZZ s,t;

//Соберём невыродившиеся строки с интересующими нас показателями в отдельную матрицу
//Эта матрица будет находится на процессе 0
//И первыми двумя её строками будут показатели при разложениях гладких а и b

//Подсчитаем количество неисключённых и невыродившихся строк в нашем стрипе матрицы
bool will_count;
int good_lines=0;

for (i=0;i<base_ub;i++)
  if (!processed_vectors[i]) {
    will_count=false;
    for (j=0;j<d1_size;j++)
	if (d1_must_leave[j]&&(base[i][j]!=0))
	  will_count=true;
    if (will_count) 
	good_lines++;
    }

int all_good_lines=0;
MPI_Allreduce(&good_lines,&all_good_lines,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);


if (all_good_lines==0) {
	if (crank==0) {
	  cout << "В матрице нет невыродившихся строк" << endl;
	  cout << "Рекомендуется увеличить размер базиса" << endl;
	  }
	MPI_Finalize();
	exit(0);
  }


ZZ ** new_base;
long new_base_ub=0;
if (crank==0) {
  cout << "all_good_lines = " << all_good_lines << endl;
  cout << "d1_must_leave_nonzero = " << d1_must_leave_nonzero << endl;

  //Создадим матрицу нужного размера
  new_base = new ZZ*[all_good_lines+2];
  for (i=0;i<all_good_lines+2;i++) {
      new_base[i]=new ZZ[d1_must_leave_nonzero];
      }
  //Запишем в неё разложения a и b
    j=0;
  for (i=0;i<d1_size;i++) {
    if (d1_must_leave[i]) {
      new_base[0][j]=a_smooth_decomp[i];
      new_base[1][j]=b_smooth_decomp[i];
      j++;
      }
    }
  new_base_ub=2;

  //Запишем строки, принадлежащие нашему процессу
  for (i=0;i<base_ub;i++)
    if (processed_vectors[i]==false) {
      will_count=false;
      for (j=0;j<d1_size;j++)
	  if (d1_must_leave[j]&&(base[i][j]!=0))
	    will_count=true;
      if (will_count) {
	  k=0;
	  for (j=0;j<d1_size;j++) {
	    if (d1_must_leave[j]) {
	      new_base[new_base_ub][k]=base[i][j];
	      k++;
	      }
	    }
	  new_base_ub++;
	  }
      }

  //Примем строки, принадлежащие другим процессам
  for (k=1;k<csize;k++) {
    int cur_count;
    MPI_Status status;
    MPI_Recv(&cur_count,1,MPI_INT,k,888,MPI_COMM_WORLD,&status);
    for (i=0;i<cur_count;i++) {
      for (j=0;j<d1_must_leave_nonzero;j++) {
	MPI_NTL_Recv(new_base[new_base_ub][j],k);
	}
      new_base_ub++;
      }
    }
  }
else {
  //Отошлём строки главному процессу
  MPI_Send(&good_lines,1,MPI_INT,0,888,MPI_COMM_WORLD);

  for (i=0;i<base_ub;i++)
    if (processed_vectors[i]==false) {
      will_count=false;
      for (j=0;j<d1_size;j++)
	  if (d1_must_leave[j]&&(base[i][j]!=0))
	    will_count=true;
      if (will_count) {
	  for (j=0;j<d1_size;j++) 
	    if (d1_must_leave[j]) {
	      MPI_NTL_Send(base[i][j],0);
	      }
	  }
      }
  }

MPI_Barrier(MPI_COMM_WORLD);
ZZ alpha,beta;
if (crank==0) {

  bool result;
  result = solve_matrix(new_base, all_good_lines+2, d1_must_leave_nonzero, r, //in
	    alpha, beta); //out
  if (!result) {
    cout << "Функция solve_matrix не смогла найти нетривиальные коэффициенты alpha, beta" << endl;
    }
  else {
    cout << "alpha = " << alpha << " beta = " << beta << endl;
    }
  }

s=(a_smooth_ind*alpha)%r;
t=(b_smooth_ind*beta)%r;

//Шаг 3 - решение показательного уравнения

if (crank==0) {
	ZZ divisor,t_norm,ans,corr;
	ZZ a_power=s;
	ZZ b_power=t;
	cout << "a_power = " << a_power << "; " << "b_power = " << b_power << endl;

	if (a_power==0 && b_power==0) {
		cout << "Сравнение выродилось, решение невозможно" << endl;
		cout << "Рекомендуется увеличить размер базиса" << endl;
		MPI_Finalize();
		exit(0);
		}
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
MPI_Finalize();
return 0;
}
