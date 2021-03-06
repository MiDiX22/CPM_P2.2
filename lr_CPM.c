#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define frand(M) (M*(((double)rand())/RAND_MAX))

#define N 5000000  

double X[N];
double Y[N];
/* rank del proces    */
int    el_meu_rank;
/* numero de processos        */
int    world_size;

double cost (int nn, double vx[], double vy[], double t0, double t1)
 {
  int i;
  double val,sum=0.0, sum_global;

  for(i=(el_meu_rank)*(nn/world_size);i<(nn/world_size)*(el_meu_rank+1);i++)
  {
    val = t0 + t1*vx[i] - vy[i];
    sum += val * val;
  }
  MPI_Reduce(&sum, &sum_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  sum_global /= 2*nn;
  return(sum_global);
}

int gradientDescent (int nn, double vx[], double vy[], double alpha, double *the0, double *the1)
 {
  int i;
  double val;
  double z0, z1, z0_global, z1_global;
  double c=0,ca;
  double t0=*the0, t1=*the1;
  double a_n = alpha/nn;
  int iter = 0, condition;
  double error = 0.000009; // cinc decimals
  
  do
  {
    z0 = z1 = 0.0;
    for(i=(el_meu_rank)*(nn/world_size);i<(nn/world_size)*(el_meu_rank+1);i++)
     {
      val = t0 + t1*vx[i] - vy[i];
      z0 += val;
      z1 += val * vx[i];
     }
    MPI_Allreduce(&z0, &z0_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&z1, &z1_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t0 -= z0_global * a_n;
    t1 -= z1_global * a_n;
    iter++;
    ca = c;
    c = cost(nn,vx,vy,t0,t1);
    condition = fabs(c - ca) > error;
    MPI_Bcast(&condition,1,MPI_INT,0,MPI_COMM_WORLD);
  }
  while (condition);
  *the0 = t0;
  *the1 = t1;
  return(iter);
}

int main()
{
  /* Inicialitzar MPI */
  MPI_Init(NULL, NULL);
  /* Obtenir el rank del proces  */
  MPI_Comm_rank(MPI_COMM_WORLD, &el_meu_rank);
  /* Obtenir el numero total de processos */
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int i;
  double ct;
  double theta0=0, theta1=1;

  srand(1);
  for (i=0;i<N;i++) 
   {
    X[i] = frand(13);
    Y[i] = frand(9) + ((1.66 + (frand(0.9))) *  X[i]) * X[i] ;
   }

  //for (i=0;i<N;i++) printf("%g %g\n",X[i],Y[i]);
 
  i=gradientDescent (N, X, Y, 0.01, &theta0, &theta1);
  ct=cost(N,X,Y,theta0,theta1);
  if (el_meu_rank == 0)
  {
    printf ("(%d) Theta; %g, %g  cost: %g\n",i,theta0,theta1,ct);
  }

  MPI_Finalize();
  return(0);
}
