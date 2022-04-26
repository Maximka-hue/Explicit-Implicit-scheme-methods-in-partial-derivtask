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

int main(void/*int argc, char * argv[]*/)
{
    constexpr static double tt{static_cast<double>(X/dsteps)};
    constexpr static double hh{static_cast<double>(Y/dsteps)};
    ofstream config; 
    ofstream file_add; 
    config.open("ConfigLayers.txt",ios::trunc);
    file_add.open("Results3_21.txt",ios::trunc);
    using V = vector<vector<double>>;
    using Vv = valarray<double>;
    vector<Vv> Vec_vals;
    Vec_vals.reserve(dsteps);
    for(int i=0;i<dsteps;i++){
    Vec_vals[i].resize(dsteps);
    Vec_vals[i] = 0;
    }
   // file_add<<"This nice, vector with all nulls,last is \t"<<Vec_vals[2][dsteps-1]<<"\n";
///Сейчас все вектора нули размера dsteps*dsteps
    int new_layer = 2;
    Vv vu;
    vu.resize(dsteps);
    vu=0;
    Vv vc;
    vc.resize(dsteps);
    vc=0;
    Vv vb;
    vb.resize(dsteps);
    vb=0;
    for(;new_layer<dsteps;new_layer++){
    vu[0]= tt * sin(new_layer*tt) + tt * (vc[1] - vc[0]) / hh + vc[0];
    file_add<<" first\t "<<vu[0]<<"\n";// Вычислили первый элемент на верхнем слое
    for(int new_step = 1;new_step<dsteps-1;new_step++){
   // double wi = 1 - 16 * hh * hh * (1 - hh) * (1 - hh);//x=1
   // vu[1] = (vc[1] * (hh / tt + 2 * tt / hh) + 
     //            vc[0] * ((hh * hh) / (tt * tt) - hh / tt - 2 * tt / hh + 2 * wi * tt * tt) + 
       //                      hh * hh / tt * sin(1 * tt) + 2 * tt * sin(1 * tt))+vb[0] ; //x=1, y = 2nd cell
    double x = new_step*hh;
    double w = 1 - 16 * x * x * (1 - x) * (1 - x);
    if(new_layer ==3&&new_step ==1){file_add<<vc[new_step - 1]<<"_________________/n";}
    vu[new_step] = tt * tt / (hh * hh) * (vc[new_step + 1] - 2 * vc[new_step] 
                            + vc[new_step - 1]) - 2 * w * tt * tt * vc[new_step] + 2 * vc[new_step] - vb[new_step];
    if(new_step%25 ==0){file_add<<vu[new_step]<<"\tat step"<<new_step<<"\n";}
    //vu[dsteps-2] = vc[dsteps-2] + tt / hh * (vc[dsteps-1] - vc[dsteps-2]);
    if(new_step==dsteps-2){vu[dsteps-1]=vc[dsteps-1]+tt * (vc[dsteps-2] - vc[dsteps-1]) / hh;
    file_add<<"Last horizontal step"<<vu[dsteps-1]<<"\n";}//закончили слой по оси x
    //if(new_layer==4 && new_step==2){file_add<<vu[3]<<"  forth layer third element \n";}
            }
    file_add<<"Finish layer #"<<new_layer<<"\n";
    //if(new_layer==dsteps-1){for(int i=0;i<dsteps;i++){file_add<<"  "<<vu[i];if(i%8==0){file_add<<"\n";};}}
    if(new_layer==dsteps-1){file_add<<"This is exact result";
        for(int i=0;i<dsteps;i++){
        file_add<<vu[i]<<" ";
        if(i%10 ==0){file_add<<"\n";}
        }    Vec_vals[0]=vu;
    }
    swap_ranges(std::begin(vu),end(vu),std::begin(vc));
    swap_ranges(std::begin(vu),end(vu),std::begin(vb));
    vu = 0;
    //config<<"This is bottom layer\t"<<vb[0]<<"Changing members up and center"<<vc[0]<<"/n This is up"<<vu[0]<<"\n";
    }
    //_______________________________________________________________________________________________________________
    double nnew_layer = 2;
    Vv vuu;
    vuu.resize(ddsteps);
    vuu=0;
    Vv vcc;
    vcc.resize(ddsteps);
    vcc=0;
    Vv vbb;
    vbb.resize(ddsteps);
    vbb=0;
    for(;nnew_layer<ddsteps;nnew_layer++){
    vuu[0]= tt * sin(nnew_layer*tt) + tt * (vcc[1] - vcc[0]) / hh + vcc[0];
    config<<" first\t "<<vuu[0]<<"\n";// Вычислили первый элемент на верхнем слое
    for(int nnew_step = 1;nnew_step<ddsteps-1;nnew_step++){
   // double wi = 1 - 16 * hh * hh * (1 - hh) * (1 - hh);//x=1
   // vu[1] = (vc[1] * (hh / tt + 2 * tt / hh) + 
     //            vc[0] * ((hh * hh) / (tt * tt) - hh / tt - 2 * tt / hh + 2 * wi * tt * tt) + 
       //                      hh * hh / tt * sin(1 * tt) + 2 * tt * sin(1 * tt))+vb[0] ; //x=1, y = 2nd cell
    double xx = nnew_step*hh;
    double ww = 1 - 16 * xx * xx * (1 - xx) * (1 - xx);
    //if(new_layer ==3&&new_step ==1){file_add<<vc[new_step - 1]<<"_________________/n";}
    vuu[nnew_step] = tt * tt / (hh * hh) * (vcc[nnew_step + 1] - 2 * vcc[nnew_step] 
                            + vcc[nnew_step - 1]) - 2 * ww * tt * tt * vcc[nnew_step] + 2 * vcc[nnew_step] - vbb[nnew_step];
    if(nnew_step%25 ==0){config<<vuu[nnew_step]<<"\tat step"<<nnew_step<<"\n";}
    //vu[dsteps-2] = vc[dsteps-2] + tt / hh * (vc[dsteps-1] - vc[dsteps-2]);
    if(nnew_step==ddsteps-2){vuu[ddsteps-1]=vcc[ddsteps-1]+tt * (vcc[ddsteps-2] - vcc[ddsteps-1]) / hh;
    config<<"Last horizontal step"<<vuu[ddsteps-1]<<"\n";}//закончили слой по оси x
    //if(new_layer==4 && new_step==2){file_add<<vu[3]<<"  forth layer third element \n";}
            }
    config<<"Finish layer #"<<nnew_layer<<"\n";
    //if(new_layer==dsteps-1){for(int i=0;i<dsteps;i++){file_add<<"  "<<vu[i];if(i%8==0){file_add<<"\n";};}}
    if(nnew_layer==ddsteps-1){config<<"This is exact result";
        for(int i=0;i<ddsteps;i++){
        config<<vuu[i]<<" ";
        if(i%10 ==0){config<<"\n";}
        }    Vec_vals[1]=vuu;
    }
    swap_ranges(std::begin(vuu),end(vuu),std::begin(vcc));
    swap_ranges(std::begin(vuu),end(vuu),std::begin(vbb));
    vuu = 0;
    //config<<"This is bottom layer\t"<<vb[0]<<"Changing members up and center"<<vc[0]<<"/n This is up"<<vu[0]<<"\n";
    }
    //Нормы
    double max = Vec_vals[0].max();
    double Max = Vec_vals[1].max();
    double Delta = abs(Vec_vals[1]-Vec_vals[0]).max();
    double delta = Delta/max;
    file_add<<"\n"<<dsteps<<" & "<<ddsteps<<"\n";
    file_add<<"\n Максимум погрешности C_h  "<<Delta<<"\n";
    double Delta2 =Max-max;
    file_add<<"\nПогрешность максимумов "<<Delta2<<"\n"<<"delta  "<<delta<<"\n";
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
    file_add<<"\nНорма L_h  "<<Norma<<"\n"<<"delta  "<<Reference<<"\n";

    file_add.close();
    config.close();}
   

/*
200 & 250
Y=10*x;
 Максимум погрешности C_h  0.104669

Погрешность максимумов 0.104669
delta  0.459699

Норма L_h  0.0695623
delta  0.825868
80 & 250

 Максимум погрешности C_h  6639.33

Погрешность максимумов 6639.11
delta  29903.2

Норма L_h  3289.21
delta  42699.8
100 * 250

 Максимум погрешности C_h  1.03194

Погрешность максимумов 0.809948
delta  4.64864

Норма L_h  1.28284
delta  16.7129

150 & 250

 Максимум погрешности C_h  0.28095

Погрешность максимумов 0.28095
delta  1.26603

Норма L_h  0.227179
delta  2.9742

100 & 200

 Максимум погрешности C_h  0.407271

Погрешность максимумов 0.407271
delta  1.83467

Норма L_h  0.403321
delta  5.2545

100 & 300

 Максимум погрешности C_h  692.668

Погрешность максимумов 692.446
delta  3120.31

Норма L_h  359.532
delta  4684.01

200 & 250

 Максимум погрешности C_h  0.103005

Погрешность максимумов 0.103005
delta  0.464262

Норма L_h  0.0656159
delta  0.861194
*/