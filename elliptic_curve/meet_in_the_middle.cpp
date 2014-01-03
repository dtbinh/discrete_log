/*
Распределённая программа для вычисления дискретного логарифма на эллиптической кривой
методом встречи посередине
Автор: Игорь Сидоров, аспирант каф. БИТ ТТИ ЮФУ
*/
#include "ec_point.h"
#include "ec_routine.h"
#include "mpi_ec_routines.h"
#include <time.h>
#include <NTL/ZZ.h>
#include <iostream>
#include <fstream>
#include "mpi.h"
using namespace std;

long ** round_robin_scheme;
long rounds;

int main(int argc, char ** argv) {
//Инициализируем MPI
MPI_Init(&argc,&argv);

//Узнаём кто мы
int crank,csize;
MPI_Comm_rank(MPI_COMM_WORLD,&crank);
MPI_Comm_size(MPI_COMM_WORLD,&csize);

SetSeed(to_ZZ(time(NULL))+crank*1000);

//Инициализируем схему обмена между процессами
init_round_robin_scheme(csize,round_robin_scheme,rounds);

long i,j,k;

//Чтение исходных данных
//Ввести данные
//Кривую
ifstream ftask("task.txt");
ZZ A,B;
ftask >> A >> B;
//Модуль
ZZ p;
ftask >> p;
ec_point::Init(A,B,p);

//Образующую точку
ZZ x,y;
ftask >> x >> y;
ec_point Q= to_ec_point(x,y);
//Порядок точки
ZZ r;
ftask >> r;
//Сама точка
ftask >> x >> y;
ec_point P=to_ec_point(x,y);
ftask.close();

//Считываем параметры алгоритма
ifstream fparam("param-mitm.txt");
long sum_db_size;
long block_size;
long add_strip_size;

//Размер базы данных
fparam >> sum_db_size;
fparam >> block_size;
fparam >> add_strip_size;
fparam.close();

double time_2,time_1;

time_1=MPI_Wtime();

//Задача P=Q*l
//Найти l (mod r)

//Создаём буфер
ec_point * temp_buf = new ec_point[block_size];

//Создадим буфер для пересылок
int res_size = block_size*5*(NumBytes(p)+1);
unsigned char * bin_buf = new unsigned char [res_size];
int bin_buf_ub,bin_buf_size;

//Массив курсоров
long * temp_buf_cursors = new long[block_size];

//Счётчик
long * temp_buf_counter = new long[csize];
long * temp_buf_last_cursor = new long[csize];

//База данных этого процесса

long max_db_size = sum_db_size/csize+add_strip_size;
ec_point * data_base = new ec_point[max_db_size];
long data_base_ub=0;
long current_sum_db_size;

do {

  //Обнулим буфера
  for (i=0;i<csize;i++) {
	  temp_buf_counter[i]=0;
	  temp_buf_last_cursor[i]=-1;
	  }

  //Пока не наберём нужное число точек в буфере
  for (i=0;i<block_size;i++) {
	  //Генерируем точку
	  ZZ lj= RandomBnd(r-1)+1;
  
	  ec_point ljQ = Q*lj;
	  //Устанавливаем показатели
	  ljQ.auxQ = lj;
	  ljQ.auxP = to_ZZ(0);
	  //помещаем точку в буфер
	  temp_buf[i]=ljQ;
	  //Указываем, какому процессу должна принадлежать данная точка
	  long process_id=to_long(ljQ.get_x()/(p/to_ZZ(csize)));
	  //Подправим, если мы ненароком вышли за границы массива
	  if (process_id>=csize) 
		  process_id=csize-1;
	  
	  //Инкрементим счётчик
	  temp_buf_counter[process_id]++;
	  //Устанавливаем курсоры
	  temp_buf_cursors[i]=temp_buf_last_cursor[process_id];
	  temp_buf_last_cursor[process_id]=i;
	  }

  //Передадим элементы сами себе
  for (j=temp_buf_last_cursor[crank];j>-1;j=temp_buf_cursors[j]) {
    if (data_base_ub<max_db_size) {
      data_base[data_base_ub++]=temp_buf[j];
      }
    }

  //Теперь начнём рассылку
  for (i=0;i<rounds;i++) {
    long partner_id = round_robin_scheme[crank][i];
    if (partner_id==-1)
      continue;
    if (crank>partner_id) {
      //Сначала передача, потом приём
      //Отошлём размер передаваемых данных
      //Затем сами данные
      MPI_Send(&(temp_buf_counter[partner_id]),1,MPI_LONG,partner_id,99,MPI_COMM_WORLD);
      //Упакуем в буфер
      bin_buf_ub=0;
      for (j=temp_buf_last_cursor[partner_id];j>-1;j=temp_buf_cursors[j]) {
	EC_Pack(temp_buf[j],bin_buf,bin_buf_ub);
	}
      Buf_send(bin_buf,bin_buf_ub,partner_id);
      //Примем размер передаваемых данных
      //Затем сами данные
      long received_len;
      MPI_Status mpi_status;
      MPI_Recv(&received_len,1,MPI_LONG,partner_id,99,MPI_COMM_WORLD,&mpi_status);
      ec_point t_point;
      Buf_recv(bin_buf, bin_buf_size, partner_id);
      bin_buf_ub=0;
      for (j=0;j<received_len;j++) {
	  EC_Unpack(t_point,bin_buf,bin_buf_ub,bin_buf_size);
	  if (data_base_ub<max_db_size) {
	    data_base[data_base_ub++]=t_point;
	    }
	  }
      }
    else {
      //Сначала приём, затем передача
      //Примем размер передаваемых данных
      //Затем сами данные
      long received_len;
      MPI_Status mpi_status;
      MPI_Recv(&received_len,1,MPI_LONG,partner_id,99,MPI_COMM_WORLD,&mpi_status);
      ec_point t_point;
      Buf_recv(bin_buf, bin_buf_size, partner_id);
      bin_buf_ub=0;
      for (j=0;j<received_len;j++) {
	  EC_Unpack(t_point,bin_buf,bin_buf_ub,bin_buf_size);
	  if (data_base_ub<max_db_size) {
	    data_base[data_base_ub++]=t_point;
	    }
	  }
      //Отошлём размер передаваемых данных
      //Затем сами данные
      MPI_Send(&(temp_buf_counter[partner_id]),1,MPI_LONG,partner_id,99,MPI_COMM_WORLD);
      bin_buf_ub=0;
      for (j=temp_buf_last_cursor[partner_id];j>-1;j=temp_buf_cursors[j]) {
	  EC_Pack(temp_buf[j],bin_buf,bin_buf_ub);
	  }
      Buf_send(bin_buf,bin_buf_ub,partner_id);
      }
    }

  //Подсчитаем, сколько всего элементов в БД
  MPI_Allreduce(&data_base_ub,&current_sum_db_size,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  } while (current_sum_db_size<sum_db_size);


//Этап построения базы завершён
time_2=MPI_Wtime();
if (crank==0) 
  cout << "База точек размером " << current_sum_db_size << " построена за " << time_2-time_1 << " секунд" << endl;

//Этап сортировки базы
time_1=MPI_Wtime();
qsort(data_base,data_base_ub,sizeof(ec_point),ec_compare);
MPI_Barrier(MPI_COMM_WORLD);

time_2=MPI_Wtime();
if (crank==0)
  cout << "База точек отсортирована за " << time_2-time_1 << " секунд" << endl;


time_1=MPI_Wtime();
//Этап поиска в БД
long found=0,all_found=0;
ec_point * search_buf = new ec_point[block_size*2];
long search_buf_size=block_size*2;
long search_buf_ub;
ec_point point_Q,point_P;

do {
  //Заполняем буфер случайными степенями P
  //Обнулим буфера
  for (i=0;i<csize;i++) {
	  temp_buf_counter[i]=0;
	  temp_buf_last_cursor[i]=-1;
	  }
  //Пока не наберём нужное число точек в буфере
  for (i=0;i<block_size;i++) {
	  //Генерируем точку
	  ZZ lj= RandomBnd(r-1)+1;

	  ec_point ljP = P*lj;
	  //Устанавливаем показатели
	  ljP.auxQ = to_ZZ(0);
	  ljP.auxP = lj;
	  //помещаем точку в буфер
	  temp_buf[i]=ljP;
	  //Указываем, какому процессу должна принадлежать данная точка
	  long process_id=to_long(ljP.get_x()/(p/to_ZZ(csize)));
	  //Подправим, если мы ненароком вышли за границы массива
	  if (process_id>=csize) 
		  process_id=csize-1;
	  
	  //Инкрементим счётчик
	  temp_buf_counter[process_id]++;
	  //Устанавливаем курсоры
	  temp_buf_cursors[i]=temp_buf_last_cursor[process_id];
	  temp_buf_last_cursor[process_id]=i;
	  }

  //Распределяем буфер между процессами
  search_buf_ub=0;
  //Передадим элементы сами себе
  for (j=temp_buf_last_cursor[crank];j>-1;j=temp_buf_cursors[j]) {
    if (search_buf_ub<search_buf_size) {
      search_buf[search_buf_ub++]=temp_buf[j];
      }
    }

  //Теперь начнём рассылку
  for (i=0;i<rounds;i++) {
    long partner_id = round_robin_scheme[crank][i];
    if (partner_id==-1)
      continue;
    if (crank>partner_id) {
      //Сначала передача, потом приём
      //Отошлём размер передаваемых данных
      //Затем сами данные

      MPI_Send(&(temp_buf_counter[partner_id]),1,MPI_LONG,partner_id,99,MPI_COMM_WORLD);
      bin_buf_ub=0;
      for (j=temp_buf_last_cursor[partner_id];j>-1;j=temp_buf_cursors[j]) {
	  EC_Pack(temp_buf[j],bin_buf,bin_buf_ub);
	  }
      Buf_send(bin_buf,bin_buf_ub,partner_id);
      //Примем размер передаваемых данных
      //Затем сами данные
      long received_len;
      MPI_Status mpi_status;
      MPI_Recv(&received_len,1,MPI_LONG,partner_id,99,MPI_COMM_WORLD,&mpi_status);
      ec_point t_point;
      Buf_recv(bin_buf, bin_buf_size, partner_id);
      bin_buf_ub=0;
      for (j=0;j<received_len;j++) {
	  EC_Unpack(t_point,bin_buf,bin_buf_ub,bin_buf_size);
	  if (search_buf_ub<search_buf_size) {
	    search_buf[search_buf_ub++]=temp_buf[j];
	    }
	  }
      }
    else {
      //Сначала приём, затем передача
      //Примем размер передаваемых данных
      //Затем сами данные
      long received_len;
      MPI_Status mpi_status;
      MPI_Recv(&received_len,1,MPI_LONG,partner_id,99,MPI_COMM_WORLD,&mpi_status);
      ec_point t_point;
      Buf_recv(bin_buf, bin_buf_size, partner_id);
      bin_buf_ub=0;
      for (j=0;j<received_len;j++) {
	  EC_Unpack(t_point,bin_buf,bin_buf_ub,bin_buf_size);
	  if (search_buf_ub<search_buf_size) {
	    search_buf[search_buf_ub++]=temp_buf[j];
	    }
	  }
      //Отошлём размер передаваемых данных
      //Затем сами данные
      MPI_Send(&(temp_buf_counter[partner_id]),1,MPI_LONG,partner_id,99,MPI_COMM_WORLD);
      bin_buf_ub=0;
      for (j=temp_buf_last_cursor[partner_id];j>-1;j=temp_buf_cursors[j]) {
	  EC_Pack(temp_buf[j],bin_buf,bin_buf_ub);
	  }
      Buf_send(bin_buf,bin_buf_ub,partner_id);
      }
    }

  //Процессы производят поиск и обмениваются результатами поиска
  found=0;
  for (i=0;i<search_buf_ub;i++) {
    ec_point * search_result = (ec_point *) bsearch(&(search_buf[i]),data_base,data_base_ub,sizeof(ec_point),ec_compare_bs);
    if (search_result!=NULL) {
      point_Q=*search_result;
      point_P=search_buf[i];
      found=1;
      break;
      }
    }
  MPI_Allreduce(&found,&all_found,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  } while (all_found==0);

time_2=MPI_Wtime();
if (crank==0)
  cout << "Нашли пару за " << time_2-time_1 << " секунд" << endl;

if (found==1) {
    cout << "Нашли пару!" << endl;
    cout << "точка Q*" << point_Q.auxQ << "=";
    point_Q.print(cout);
    cout << "точка P*" << point_P.auxP << "=";
    point_P.print(cout);
    ZZ result = calc_answer(point_Q, point_P,Q,P,r);
    cout << "Ответ = " << result << endl;
    if (Q*result==P) {
      cout << "Ответ верен" << endl;
      }
    else {
      cout << "Ответ НЕВЕРЕН" << endl;
      }
    }

MPI_Barrier(MPI_COMM_WORLD);

MPI_Finalize();

return 0;
}
