#include <iostream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <fstream>
using namespace std;
// rectangal method
void show(double *MVal,int r){
for(int l=0; l<r; l++){
printf("%6f", MVal[l]);
}
printf("/n");
}

double integer1(double (*f)(double), double a, double b, int n=1000){
double sum=0;
int k=0;
double h=(b-a)/n;
while(k<=n){
sum+=f(a+k*h)*h;
k++;
}
return sum;
}
//iterator method
//double Irecur(double integer2,double n);
//I(n)=1/n-I(n-1);


double integer2(int iter=30,double In=log(7/6)){
for(int k=1;k<=iter;k++){
In=log(6/7);
In=1/k-In;
}
return integer2();
}

double integer3(int iter=30,double In=0){
for(int i=60;i!=iter+1;i--){
In=1/(6*i)-In/6;
}
return integer3();
}

double F(double n=1){
double x;
return x.pow(n)/(x*x*x*x*x*x+1);
}

int main(){
system("chcp 1251>nul");
cout<<"���������� ���������� p.1\n";
cout<<integer1(F, 0.0, 1.0)<<endl;
cout<<"���������� ���������� p.2\n";
cout<<integer2()<<endl;
cout<<"���������� ���������� p.2\n";
cout<<integer3()<<endl;


system("pause.nul");
return 0;
}

