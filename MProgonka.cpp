#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

int main(void){
	
   int     i, j, n, m, k;
   double  h, L, x, a, T, tau, sum=0, max=0,sum_2=0,max_2=0,n_2,m_2,h_2,tau_2,cef;
   x   = 0;
   a   = 0.01;
   L   = 1;
   T   = 1;
   n   = 50;
   m   = 50;
   h   = L/(n-1);
   tau = T/(m-1);
   
   n_2   = n;
   m_2   = 1000*m;
   h_2   = L/(n_2-1);
   tau_2 = T/(m_2-1);
   //printf ("%f", tau);
   double *u, *u_2;
   double *u_befor,*w, *A, *B, *C, *D;
   double *u_befor_befor;
   double *alpha, *beta, *gamma;
   
   
   if (!(u = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
       
   if (!(u_2 = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
    if (!(u_befor = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
    if (!(u_befor_befor = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
    if (!(w = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
      if (!(A = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
      if (!(B = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
      if (!(C = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
       if (!(D = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
     if (!(alpha = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
     if (!(beta = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
     if (!(gamma = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
    u_befor[0]=sin(tau);
    u_befor_befor[0]=0;
   
    for(j=0;j<n;j++){
    x=j*h;
    u[j]=0;
    //printf ("%f", u[j]);
	u_befor[j]=0;
    u_befor_befor[j]=0;
    w[j]=0;
    
    A[j]=0;
	B[j]=0;
	C[j]=0;
	D[j]=0;
	gamma[j]=0;
	alpha[j]=0;
	beta[j]=0;
	}
	
	for (j=2;j<m;j++){
		u[0]=sin(tau*j);
		u[n-1]=0;
		w[0]=2*a*pow((u_befor[1]-u_befor[0])/h,2);
    for (i=1;i<n-1;i++){
    	x=i*h;
        w[i]=4*x*(1-x)+2*a*pow((u_befor[i]-u_befor[i-1])/h,2);
    	
    	C[i]=-(pow(tau,2)/(pow(h,2)))*(w[i-1]);
    	B[i]=1+(pow(tau,2)/(pow(h,2)))*((w[i]+w[i-1]));
    	A[i]=-(pow(tau,2)/(pow(h,2)))*w[i];
    	D[i]=2*u_befor[i]-u_befor_befor[i];
    
    	//u[i]=(pow(tau,2)/(pow(h,2)))*((w[i]-w[i-1])*(u_befor[i+1]-u_befor[i])+w[i]*(u_befor[i+1]-2*u_befor[i]+u_befor[i-1]))+2*u_befor[i]-u_befor_befor[i];
    	//w[i]=4*x*(1-x)+2*a*pow((u_befor[i]-u_befor[i-1])/h,2);
       
    }
    
    A[0]=0;
    B[0]=1;
    C[0]=0;
    D[0]=sin(tau*j);
    
    alpha[0]=0;
    beta[0]=sin(tau*j);
    
    alpha[1]=-C[0]/B[1];
    beta[1]=D[1]/B[1];
   
    
    	for (i=0;i<n-1;i++){
    	gamma[i] =1/(B[i+1]+C[i+1]*alpha[i]);
    	alpha[i+1] = -A[i+1]*gamma[i];
        beta[i+1]  = gamma[i]*(D[i+1]-C[i+1]*beta[i]); 
         }
         
 
        for (i=n-2;i>0;i--){
        u[i]=alpha[i]*u[i+1]+beta[i];
        }
        
         for (i=0;i<n;i++){
    	u_befor_befor[i]=u_befor[i];
       	u_befor[i]=u[i];
    }
         
    
    //printf ("\n\n\n");

}
 for (i=0;i<n;i++){
    u_2[i]=u[i];
   }
    

   n   = n_2;
   m   = m_2;
   h   = h_2;
   tau = tau_2;

u_befor[0]=sin(tau);
    u_befor_befor[0]=0;
   
    for(j=0;j<n;j++){
    x=j*h;
    u[j]=0;
    //printf ("%f", u[j]);
	u_befor[j]=0;
    u_befor_befor[j]=0;
    w[j]=0;
    
    A[j]=0;
	B[j]=0;
	C[j]=0;
	D[j]=0;
	gamma[j]=0;
	alpha[j]=0;
	beta[j]=0;
	}

	for (j=2;j<m;j++){
		u[0]=sin(tau*j);
		u[n-1]=0;
		w[0]=2*a*pow((u_befor[1]-u_befor[0])/h,2);
    for (i=1;i<n-1;i++){
    	
    	
    	x=i*h;
        w[i]=4*x*(1-x)+2*a*pow((u_befor[i]-u_befor[i-1])/h,2);
    	
    	C[i]=-(pow(tau,2)/(pow(h,2)))*(w[i-1]);
    	B[i]=1+(pow(tau,2)/(pow(h,2)))*((w[i]+w[i-1]));
    	A[i]=-(pow(tau,2)/(pow(h,2)))*w[i];
    	D[i]=2*u_befor[i]-u_befor_befor[i];
    	
    
    	//u[i]=(pow(tau,2)/(pow(h,2)))*((w[i]-w[i-1])*(u_befor[i+1]-u_befor[i])+w[i]*(u_befor[i+1]-2*u_befor[i]+u_befor[i-1]))+2*u_befor[i]-u_befor_befor[i];
    	//w[i]=4*x*(1-x)+2*a*pow((u_befor[i]-u_befor[i-1])/h,2);
   

       
       
    }
    
    A[0]=0;
    B[0]=1;
    C[0]=0;
    D[0]=sin(tau*j);
    
    alpha[0]=0;
    beta[0]=sin(tau*j);
    
    alpha[1]=-C[0]/B[1];
    beta[1]=D[1]/B[1];
   
    
    	for (i=0;i<n-1;i++){
    	gamma[i] =1/(B[i+1]+C[i+1]*alpha[i]);
    	alpha[i+1] = -A[i+1]*gamma[i];
        beta[i+1]  = gamma[i]*(D[i+1]-C[i+1]*beta[i]); 
         }
         
   
         
        for (i=n-2;i>0;i--){
        u[i]=alpha[i]*u[i+1]+beta[i];
        }
        
         for (i=0;i<n;i++){
    	u_befor_befor[i]=u_befor[i];
       	u_befor[i]=u[i];
    }
         
    
    
}

    printf ("\n\n\n");
   for (i=0;i<n;i++){
   	printf ("%f\n", u_2[i]);
   }
   
       printf ("\n\n\n");
   for (i=0;i<n;i++){
   	printf ("%f\n", u[i]);
   }
   
   
    max=fabs(u[0]-u_2[0]);
    for(i=0;i<n;i++){
   	if (fabs(u[i]-u_2[i])>=max)
   	max=fabs(u[i]-u_2[i]);
   }
   
    printf ("\n\n\n");
    printf ("%f\n", max);
    
    sum=0;
    for(i=0;i<n;i++){
    sum=sum+fabs(u[i]-u_2[i]);
    }
    
    printf ("\n\n\n");
    printf ("%f\n", sum);
    
    max_2=fabs(u_2[1]);
    for(i=0;i<n;i++){
   	if (fabs(u_2[i])>=max_2)
   	max_2=fabs(u_2[i]);
   }
   
   max=max/max_2;
   
   printf ("\n\n\n");
   printf ("%f\n", max);
   
    sum_2=0;
    for(i=0;i<n;i++){
    sum_2=sum_2+fabs(u_2[i]);
    }
    
    sum=sum/sum_2;
    
    printf ("\n\n\n");
    printf ("%f\n", sum);
   
   

  return 1;
}
