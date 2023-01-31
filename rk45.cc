#include <stdio.h>
#include <math.h>
#include "rk45.h"

const double h = 2.0e-3;
const int max_steps = 1e5;
const double gamma0 = 10.0;
const double b_imp = 50.0; //impact parameter
const double x0 = 2000.0; //initial separation in fm

int main()
{
  FILE *data = fopen("data.txt","w");

  //rk variables
  int i; //row number
  int j; //column number
  int n; //equation number
  int step_num=0; //step number
  double k[M][stage_num],alpha[M][stage_num],beta[M];

  //dynamical variables
  double tau = 0.0;
  double X[M];
  double u_squared;

  //initial conditions
  X[0] = 0.0; //t
  X[1] = x0; //x
  X[2] = b_imp; //y
  X[3] = gamma0; //gamma
  X[4] = -sqrt(gamma0*gamma0-1.0); //ux
  X[5] = 0.0; //uy

  do
    {
      //error check with u^2=1
      u_squared = fabs(X[3]*X[3] - X[4]*X[4] - X[5]*X[5] - 1.0);

      //print solution and error
      fprintf(data,"%.10f   %.10f   %.10f   %.10f   %.10f \n",X[0],X[1],X[2],X[3],u_squared);

      //calculate integration step with rk45     
      for(i=0;i<stage_num;++i)
	{
	  for(n=0;n<M;++n)
	    {
	      alpha[n][i] = 0.0;
	      for(j=0;j<i;++j)
		{
		  alpha[n][i] = alpha[n][i] + A(i,j)*k[n][j];
		}
	    }
	  for(n=0;n<M;++n)
	    {
	      k[n][i] = h*stage(n,i,X,alpha);
	    }
	}
      for(n=0;n<M;++n)
	{
	  beta[n] = 0.0;
	  for(i=0;i<stage_num;++i)
	    {
	      beta[n] = beta[n] + b(i)*k[n][i];
	    }
	  X[n] = X[n] + beta[n];
	}
      tau = tau + h;
      ++step_num;
    }
  
  while(step_num<max_steps);
  fclose(data);
  return 0;
}
		     
