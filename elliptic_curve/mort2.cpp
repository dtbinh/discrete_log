/*
Распределённая программа для вычисления дискретного логарифма на эллиптической кривой
методом встречи на случайном дереве
Автор: Игорь Сидоров, аспирант каф. БИТ ТТИ ЮФУ
*/

#include <iostream>
#include <fstream>
using namespace std;

#include "ec_routine.h"
#include "rb_tree.h"
#include "mpi_ec_routines.h"
#include "mpi.h"

rb_tree data_base;
long data_base_size;
long data_base_ub;

long ** round_robin_scheme;
long rounds;

ec_point P,Q;
ZZ one_third,two_third;
ZZ r;

void tau(ec_point& X) {
ZZ auxQ,auxP;
if (X.get_x()<one_third) {
	auxQ=X.auxQ+1;
	auxP=X.auxP;
	X=X+Q;
	X.auxQ = auxQ%r;
	X.auxP = auxP%r;
	}
else if (X.get_x()>one_third && X.get_x()<two_third) {
	auxQ=2*X.auxQ;
	auxP=2*X.auxP;
	X=X*to_ZZ(2);
	X.auxQ = auxQ%r;
	X.auxP = auxP%r;
	}
else {
	auxQ=X.auxQ+2;
	auxP=X.auxP;
	X=X+Q*to_ZZ(2);
	X.auxQ = auxQ%r;
	X.auxP = auxP%r;
	}
}

int main(int argc, char ** argv) {
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
Q= to_ec_point(x,y);
//Порядок точки
ftask >> r;
//Сама точка
ftask >> x >> y;
P=to_ec_point(x,y);
ftask.close();

//Считываем параметры алгоритма
ifstream fparam("param-mort.txt");
long sum_db_size;
long block_size;
long add_strip_size;
long k_tau;
//Размер базы данных
fparam >> sum_db_size;
fparam >> block_size;
fparam >> add_strip_size;
fparam >> k_tau;
fparam.close();

double time_2,time_1;

data_base_size = sum_db_size/csize + add_strip_size;
data_base_ub=0;

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


ec_point point_1, point_2;
long found,all_found;
long total_db_size=0;
bool hit_ceil=false;

time_1=MPI_Wtime();
do {
  found=0;
  all_found=0;
  
  //Обнулим буфера
  for (i=0;i<csize;i++) {
	  temp_buf_counter[i]=0;
	  temp_buf_last_cursor[i]=-1;
	  }

  //Пока не наберём нужное число точек в буфере
  for (i=0;i<block_size;i++) {
	  //Генерируем точку
	  ZZ lj=RandomBnd(r-1)+1;
	  ZZ li=RandomBnd(r-1)+1;
	  if (total_db_size>=sum_db_size)
	    li=0;
	  ec_point random_point = Q*lj+P*li;
	  //Устанавливаем показатели
	  random_point.auxQ = lj;
	  random_point.auxP = li;
	  //Выполним k отображений tau
	  for (j=0;j<k_tau;j++)
	    tau(random_point);
	  //помещаем точку в буфер
	  temp_buf[i]=random_point;
	  //Указываем, какому процессу должна принадлежать данная точка
	  long process_id=to_long(random_point.get_x()/(p/to_ZZ(csize)));
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
      ec_point t_point=temp_buf[j];
      Node * t_node = data_base.findNode(t_point);
      if (t_node!=NULL && (total_db_size>=sum_db_size)) {
	if (t_node->data.auxP!=t_point.auxP || t_node->data.auxQ!=t_point.auxQ) {
	  point_1=t_node->data;
	  point_2=t_point;
	  found=1;
	  }
	}
      else {
	if (data_base_ub<data_base_size) {
	  data_base.insertNode(temp_buf[j]);
	  data_base_ub++;
	  }
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
	  Node * t_node = data_base.findNode(t_point);
         if (t_node!=NULL && (total_db_size>=sum_db_size)) {
	    if (t_node->data.auxP!=t_point.auxP || t_node->data.auxQ!=t_point.auxQ) {
	      point_1=t_node->data;
	      point_2=t_point;
	      found=1;
	      }
	    }
	  else {
	    if (data_base_ub<data_base_size) {
	      data_base.insertNode(temp_buf[j]);
	      data_base_ub++;
	      }
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
	  Node * t_node = data_base.findNode(t_point);
         if (t_node!=NULL && (total_db_size>=sum_db_size)) {
	    if (t_node->data.auxP!=t_point.auxP || t_node->data.auxQ!=t_point.auxQ) {
	      point_1=t_node->data;
	      point_2=t_point;
	      found=1;
	      }
	    }
	  else {
	    if (data_base_ub<data_base_size) {
	      data_base.insertNode(temp_buf[j]);
	      data_base_ub++;
	      }
	    }
	  }
      //Отошлём размер передаваемых данных
      //Затем сами данные
      MPI_Send(&(temp_buf_counter[partner_id]),1,MPI_LONG,partner_id,99,MPI_COMM_WORLD);
      //Упакуем в буфер
      bin_buf_ub=0;
      for (j=temp_buf_last_cursor[partner_id];j>-1;j=temp_buf_cursors[j]) {
	EC_Pack(temp_buf[j],bin_buf,bin_buf_ub);
	}
      Buf_send(bin_buf,bin_buf_ub,partner_id);
      }
    }

  //Проверим, не нашёл ли кто-нибудь пару
  MPI_Allreduce(&found,&all_found,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&data_base_ub,&total_db_size,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  if (!hit_ceil && (total_db_size>=sum_db_size)) {
    time_2=MPI_Wtime();
    if (crank==0)
      cout << "Время заполнения БД: " << time_2 - time_1 << " секунд" << endl;
    time_1=MPI_Wtime();
    hit_ceil=true;
    }
  } while (all_found==0);

time_2=MPI_Wtime();
if (crank==0)
  cout << "Время поиска: " << time_2-time_1 << " секунд" << endl;


MPI_Allreduce(&data_base_ub,&total_db_size,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
if (crank==0)
    cout << "Суммарный размер БД = " << total_db_size << endl;
    
if (found==1) {
  cout << "Точка #1: Q* "<< point_1.auxQ << " + P* "<< point_1.auxP << " = ";
  point_1.print(cout);

  cout << "Точка #2: Q* "<< point_2.auxQ << " + P* "<< point_2.auxP << " = ";
  point_2.print(cout);
  ZZ result = calc_answer(point_1,point_2,Q,P,r);
  cout << "Ответ = " << result << endl;
  if (Q*result==P) {
    cout << "Ответ верен" << endl;
    }
  else {
    cout << "Ответ НЕВЕРЕН" << endl;
    }
  }

MPI_Finalize();
return 0;
}
