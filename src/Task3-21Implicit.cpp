//clang++ -std=c++17 | -std=gnu++17 /home/maximkasm/Desktop/C++/Progs/
//types
#include <array>
#include <valarray>
#include <vector>
#include <map>
#include <cstring>
#include <cwchar>
//IO
#include <iostream>
#include <iomanip>
//std
#include <cstdlib>
#include <stdexcept>//domain_error,invalid_argument
/*For open file*/
#include <cerrno>
#include <cstdio>//fopen_s,perror, /* printf */
#include <cstdint>//wchar,wint
#include <climits>//char
#include <fstream>
//time
#include <ctime>
#include <clocale>
//etc
#include <iterator>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>

#define X 1.
#define Y 1.//Domain coord
#define dsteps 150
#define ddsteps 250
using namespace std;

template<typename Matrix_type>
void LU1(vector<vector<Matrix_type>> A,vector<vector<Matrix_type>> &L,
            vector<vector<Matrix_type>> &U,int n);
void show(vector <vector <double>> A, int n);

int main(void/*int argc, char * argv[]*/)
{
    constexpr static double tt{static_cast<double>(X/dsteps)};
    constexpr static double hh{static_cast<double>(Y/dsteps)};
    ofstream inapparent; 
    inapparent.open("Implicit1.txt",ios::trunc);
    ofstream inapparentt; 
    inapparentt.open("Implicit2.txt",ios::trunc);
    //Вспомогательный вектор
    vector<vector<double>> Vec_vals;
    Vec_vals.reserve(3);
    for(int i=0;i<3;i++){
    Vec_vals[i].reserve(dsteps+1);
    fill(Vec_vals[i].begin(),Vec_vals[i].end(),0);
    //inapparent<<Vec_vals[i][dsteps-1];
    }
//Создаем матрицу А
/*
   vector<vector<double>> A;
   A.reserve(dsteps-2);
   A[0].reserve(dsteps-2);
   A[0][0] = 1 + 2*tt/hh;
   A[0][1] = -tt*tt/(hh*hh);
   fill(A[0].begin()+2,A[0].end(),0);
  // inapparent<< A[0][1]<<"KKK\n";
   for(int i=1;i<dsteps-4;i++){
    A[i].reserve(dsteps-2);
    A[i][i] = 1 + 2*tt*tt/(hh*hh);
    A[i][i-1] = -tt*tt/(hh*hh);
    A[i][i+1] = -tt*tt/(hh*hh);
    }
    inapparent<< A[dsteps-5][dsteps-4]<<"KKKaaa\n";
    A[dsteps-3][dsteps-4] =-tt*tt/(hh*hh);
    A[dsteps-3][dsteps-3] = 1 + 2*tt*tt/(hh*hh);*/
    //вектор свободных членов
    vector<double> b(dsteps,0);//n-1
    //b.reserve(dsteps);
    using Vv = valarray<double>;
    //Коэффиценты прогонки
    Vv alpha;
    alpha.resize(dsteps);
    alpha=0;
    Vv beta;
    beta.resize(dsteps);
    beta=0;
    //Векторы неизвестных
    vector<double> vb(dsteps+1,0);
    vector<double> vc(dsteps+1,0);
    vector<double> vu(dsteps+1,0);
    //v.reserve(dsteps);
    //inapparent<<b[0]<<"ll";
    //inapparent<<b[dsteps-1]<<"..."<<vb[dsteps-1]<<"\n";

    int new_layer = 2;auto ta = tt*tt/(hh*hh); auto tb = 1 + 2*tt*tt/(hh*hh);
    for(;new_layer<dsteps;new_layer++){
    //inapparent<<alpha[0] ;
   // for(;new_layer<dsteps-1;new_layer++){
    //Выражаем явно из граничных
    vu[0]= tt * sin(new_layer*tt) + tt * (vc[1] - vc[0]) / hh + vc[0];
    vu[dsteps]=vc[dsteps]+tt * (vc[dsteps-1] - vc[dsteps])/hh;
    beta[0]= vu[0]/(tb);
    alpha[0] = ta/(tb);
    for(int new_step = 1;new_step<dsteps-1;new_step++){//Уже считаем n-1 m-1
    double x = new_step*hh;
    double w = 1 - 16 * x * x * (1 - x) * (1 - x);
    //Можно было вместо коэффицентов писать соответствующие элементы матрицы,но они одинаковые,проще подставить
    b[new_step] = - 2 * w * tt * tt * vc[new_step] + 2 * vc[new_step] - vb[new_step];//Столюбец свободных
    alpha[new_step] = ta / (tb - ta* alpha[new_step - 1]);
    beta[new_step] = (b[new_step-1] + ta * beta[new_step - 1])/(tb - ta * alpha[new_step - 1]);
    //inapparent<<b[dsteps]<<beta[dsteps-2]<<beta[1]
      //          <<alpha[dsteps-2]<<alpha[dsteps-1]<<vu[dsteps-2]<<"\n";
    }
        for (int i = dsteps-2;i>0;i--){
        vu[i]=alpha[i+1]*vu[i+1]+beta[i+1];
        //inapparent<<vu[i];
        if(i%25 ==0){inapparent<<vu[i]<<"\tat step"<<i<<"\n";}
        }
    if(new_layer==dsteps-1){inapparent<<"This is exact result";
        for(int i=0;i<dsteps-1;i++){
        inapparent<<vu[i]<<" ";
        if(i%10 ==0){inapparent<<"\n";}
        }    Vec_vals[0]=vu;
    }
    swap_ranges(std::begin(vu),end(vu),std::begin(vc));
    swap_ranges(std::begin(vu),end(vu),std::begin(vb));
    fill(vu.begin(),vu.end(),0);
    }
    /*
    for(int o=0;o<dsteps-1;o++){

    for(int k = 1;k<dsteps-2;k++){
    double x = o*hh;
    double w = 1 - 16 * x * x * (1 - x) * (1 - x);
    b[k] = -2*w*
    }
    }*/
    /*
    memset(&vu,0,vu.size());
    Evalut_results.close();
    return *this;}*/
    //Вспомогательный вектор
//Создаем матрицу А
/*
   vector<vector<double>> A;
   A.reserve(dsteps-2);
   A[0].reserve(dsteps-2);
   A[0][0] = 1 + 2*tt/hh;
   A[0][1] = -tt*tt/(hh*hh);
   fill(A[0].begin()+2,A[0].end(),0);
  // inapparent<< A[0][1]<<"KKK\n";
   for(int i=1;i<dsteps-4;i++){
    A[i].reserve(dsteps-2);
    A[i][i] = 1 + 2*tt*tt/(hh*hh);
    A[i][i-1] = -tt*tt/(hh*hh);
    A[i][i+1] = -tt*tt/(hh*hh);
    }
    inapparent<< A[dsteps-5][dsteps-4]<<"KKKaaa\n";
    A[dsteps-3][dsteps-4] =-tt*tt/(hh*hh);
    A[dsteps-3][dsteps-3] = 1 + 2*tt*tt/(hh*hh);*/
    //вектор свободных членов
    vector<double> bb(ddsteps,0);//n-1
    //b.reserve(dsteps);
    //Коэффиценты прогонки
    Vv alphaa;
    alphaa.resize(ddsteps);
    alphaa=0;
    Vv betaa;
    betaa.resize(ddsteps);
    betaa=0;
    //Векторы неизвестных
    vector<double> vbb(ddsteps+1,0);
    vector<double> vcc(ddsteps+1,0);
    vector<double> vuu(ddsteps+1,0);
    //v.reserve(dsteps);
    //inapparentt<<b[0]<<"ll";
    //inapparentt<<b[dsteps-1]<<"..."<<vb[dsteps-1]<<"\n";

    new_layer = 2;
    for(;new_layer<ddsteps;new_layer++){
    //inapparentt<<alpha[0] ;
   // for(;new_layer<dsteps-1;new_layer++){
    //Выражаем явно из граничных
    vuu[0]= tt * sin(new_layer*tt) + tt * (vcc[1] - vcc[0]) / hh + vcc[0];
    vuu[dsteps]=vcc[dsteps]+tt * (vcc[dsteps-1] - vcc[dsteps])/hh;
    betaa[0]= vuu[0]/(tb);
    alphaa[0] = ta/(tb);
    for(int nnew_step = 1;nnew_step<ddsteps-1;nnew_step++){//Уже считаем n-1 m-1
    double xx = nnew_step*hh;
    double ww = 1 - 16 * xx * xx * (1 - xx) * (1 - xx);
    //Можно было вместо коэффицентов писать соответствующие элементы матрицы,но они одинаковые,проще подставить
    bb[nnew_step] = - 2 * ww * tt * tt * vcc[nnew_step] + 2 * vcc[nnew_step] - vbb[nnew_step];//Столюбец свободных
    alphaa[nnew_step] = ta / (tb - ta* alphaa[nnew_step - 1]);
    betaa[nnew_step] = (bb[nnew_step-1] + ta * betaa[nnew_step - 1])/(tb - ta * alphaa[nnew_step - 1]);
    //inapparentt<<b[dsteps]<<beta[dsteps-2]<<beta[1]
      //          <<alpha[dsteps-2]<<alpha[dsteps-1]<<vu[dsteps-2]<<"\n";
    }
        for (int i = ddsteps-2;i>0;i--){
        vuu[i]=alphaa[i+1]*vuu[i+1]+betaa[i+1];
        //inapparent<<vu[i];
        if(i%25 ==0){inapparentt<<vuu[i]<<"\tat step"<<i<<"\n";}
        }
    if(new_layer==ddsteps-1){inapparentt<<"This is exact result";
        for(int i=0;i<ddsteps;i++){
        inapparentt<<vuu[i]<<" ";
        if(i%10 ==0){inapparentt<<"\n";}
        }    Vec_vals[1]=vuu;
    }
    swap_ranges(std::begin(vuu),end(vuu),std::begin(vcc));
    swap_ranges(std::begin(vuu),end(vuu),std::begin(vbb));
    fill(vuu.begin(),vuu.end(),0);
    }

//___________________________________________//
    //vector<int>::iterator max_i;
    auto max_i = max_element(Vec_vals[0].begin(),Vec_vals[0].end());
    double max = *max_i;
    auto Max_i = max_element(Vec_vals[1].begin(),Vec_vals[1].end());
    double Max = *Max_i;
    //auto Delta_i = max_element(Vec_vals[1].begin()- Vec_vals[0].begin(),Vec_vals[1].end()- Vec_vals[0].end());
    //double Delta = *Delta_i;
    double Delta = Max - max;
    double delta = Delta/max;
    inapparent<<"\n"<<dsteps<<" & "<<ddsteps<<"\n";
    inapparent<<"\nПогрешность максимумов  C_h  "<<Delta<<"\n delta"<<delta<<"\n";
    //double Delta2 =Max-max;
    //inapparentt<<"\n Максимум погрешности"<<Delta2<<"\n"<<"delta  "<<delta<<"\n";
    double norma=0;
    for(int i =0; i<dsteps;i++){
    norma += hh*Vec_vals[0][i];
    }
    double norma2=0;
    for(int i =0; i<ddsteps;i++){
    norma2 += hh*Vec_vals[1][i];
    }
    double Norma = abs(norma - norma2);
    double Reference = Norma/norma;
    inapparent<<"\nНорма L_h  "<<Norma<<"\n"<<"delta  "<<Reference<<"\n";
    inapparent.close();
    inapparentt.close();
    return 0;
}
    
//LU Decomposition
void LU1(vector <vector <double>> A, vector <vector <double>> &L, 
		vector <vector <double>> &U, int n)
{
	U=A;

	for(int i = 0; i < n; i++)
		for(int j = i; j < n; j++)
			L[j][i]=U[j][i]/U[i][i];
	
	for(int k = 1; k < n; k++)
	{
		for(int i = k-1; i < n; i++)
			for(int j = i; j < n; j++)
				L[j][i]=U[j][i]/U[i][i];

		for(int i = k; i < n; i++)
			for(int j = k-1; j < n; j++)
				U[i][j]=U[i][j]-L[i][k-1]*U[k-1][j];
	}

}
template<typename Matrix_type>
void LU2(vector<vector<Matrix_type>> A,vector<vector<Matrix_type>> &L,
            vector<vector<Matrix_type>> &U, int n){
U=A;//upmatrix
//for(int i=0;i<n;i++){//i-raw,j-column
//static int p = 0;
for(int i =0;i<n;i++)
    transform(U[i].begin(),U[i].end(),L[i].begin(),[i](Matrix_type* x,Matrix_type* y){int p =0;x[p]/x[i];p++;});
    cout<<U[0][2]<<"Example\n";
    for(int k = 1; k < n; k++)
	{
		for(int i = k-1; i < n; i++)
			for(int j = i; j < n; j++)
				L[j][i]=U[j][i]/U[i][i];

		for(int i = k; i < n; i++)
			for(int j = k-1; j < n; j++)
				U[i][j]=U[i][j]-L[i][k-1]*U[k-1][j];
	}
}
void show(vector <vector <double>> A, int n)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			cout <<"\t"<< A[i][j] << "\t";
		}
		cout << endl;
	}
}
//Proisvedenie
void product(int M, int N, int K, const float * A, const float * B, float * C)
{
    for (int i = 0; i < M; ++i)
    {
        float * c = C + i * N;
        for (int j = 0; j < N; ++j)
            c[j] = 0;
        for (int k = 0; k < K; ++k)
        {
            const float * b = B + k * N;
            float a = A[i*K + k];
            for (int j = 0; j < N; ++j)
                c[j] += a * b[j];
        }
    }
}
//Решение уравнения By=d
template<typename Matrix_type>
vector<Matrix_type> By_d(vector<vector<Matrix_type>> B,vector<double> d,int steps){
vector<Matrix_type> y;
y.reserve(steps);
for(int j = 1; j< steps; j++){
    y.at(0) = d[0]/B[0][0];
    int sum_of_eq(0);
    //Различные способы суммирования
    for(int i(B[j].size());i>=j-1;--i)
        sum_of_eq +=B[j-1].at(i-1)*y[i-1];
    y.at(j) = (sum_of_eq+y[0])/B[j][j];
   //y.at(i) = (d.at(i) - for_each(B[i].begin(), B[i].begin()+(i-1),
     //                                           [&B,&y] (int k ){sum += k*y[k]});
}return y;}
//Решение уравнения Cx=y
template<typename Matrix_type>
void Cx_y(vector<vector<Matrix_type>> C,vector<double> y,int steps){
vector<Matrix_type> x;
x.reserve(steps);
int sum_of_c(0);
x[steps-1] = y.at(steps - 1);
for(int j = steps - 2; j>=0; --j){
    for(int i(C[j].size());i>=j-1;--i){
    sum_of_c +=C[j].at(i-1)*y[i-1];
    }
    x.at(j) = y.at(j) - (sum_of_c + x[steps-1])/C[j][j];
         //sum_of_c = accumulate(C[j].rbegin(),C[j].rbegin()+j+1, static_cast<Matrix_type>(0))
   //x[j] = y[j] - for_each(C[j].rbegin(),C[j].rbegin() + j+1, [&C,&x](int h){
  //                                                     sum_of_c += h*x[]})
}}

/*
100 & 250

Погрешность максимумов  C_h  0.156182
 delta15.0308

Норма L_h  0.0925129
delta  80.1477
100 & 250
X=10*Y
Погрешность максимумов  C_h  0
 delta-nan

Норма L_h  1.24166e+231
delta  -4.09706e+142
200 & 250

Погрешность максимумов  C_h  0.000682581
 delta0.130575

Норма L_h  0.000261452
delta  0.461834
200 & 250
X=10*Y;
Погрешность максимумов  C_h  0.000682581
 delta0.130575

Норма L_h  0.000261452
delta  0.461834
*/