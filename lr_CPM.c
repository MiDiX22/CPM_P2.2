#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define frand(M) (M*(((double)rand())/RAND_MAX))

#define N 5000000  

double X[N];
double Y[N];

double cost (int nn, double vx[], double vy[], double t0, double t1)
{
  int i;
  double val,sum=0.0;

  //bucle a paralelitzar default(none) reduction(+:sum) firstprivate(t0, t1, nn) private(i, val) shared(vx, vy) schedule(static)
  for(i=0;i<nn;i++)
    {
    val = t0 + t1*vx[i] - vy[i];
    sum += val * val;
    }
  sum /= 2*nn;
  return(sum);
}

int gradientDescent (int nn, double vx[], double vy[], double alpha, double *the0, double *the1)
{
  int i;
  double val;
  double z0,z1;
  double c=0,ca;
  double t0=*the0, t1=*the1;
  double a_n = alpha/nn;
  int iter = 0;
  double error = 0.000009; // cinc decimals

  do
    {
    z0 = z1 = 0.0;
    //Paralelitzar default(none) reduction(+:z0,z1) private(i, val) firstprivate(t0, t1, nn) shared(vx, vy) schedule(static)
    for(i=0;i<nn;i++)
      {
      val = t0 + t1*vx[i] - vy[i];
      z0 += val;
      z1 += val * vx[i];
      }
    t0 -= z0 * a_n;
    t1 -= z1 * a_n;
    iter++;
    ca = c;
    c = cost(nn,vx,vy,t0,t1);
    }
  while (fabs(c - ca) > error);
  *the0 = t0;
  *the1 = t1;
  return(iter);
}

int main()
{
  /* rank del proces    */
  int    el_meu_rank;
  /* numero de processos        */
  int    world_size;
  /* rank de l'emissor */
  int    font;
  /* rank del receptor */
  int    desti;
  /* etiqueta dels missatges */
  int    etiq = 0;
  /* espai per al missatge      */
  char missatge[100];    
  /* estat de la recepcio       */
  MPI_Status  estat;

  /* Inicialitzar MPI */
  MPI_Init(NULL, NULL);
  /* Obtenir el rank del proces  */
  MPI_Comm_rank(MPI_COMM_WORLD, &el_meu_rank);
  /* Obtenir el numero total de processos */
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (el_meu_rank == 0)
  {
    int i;
    double ct;
    double theta0=0, theta1=1;

    // sprintf(missatge, "Salutacions del proces %d ",el_meu_rank);
    // MPI_Bcast(missatge,101,MPI_CHAR,0,MPI_COMM_WORLD);

    srand(1);
    for (i=0;i<N;i++) 
      {
      X[i] = frand(13);
      Y[i] = frand(9) + ((1.66 + (frand(0.9))) *  X[i]) * X[i] ;
      }

    //for (i=0;i<N;i++) printf("%g %g\n",X[i],Y[i]);

    i=gradientDescent (N, X, Y, 0.01, &theta0, &theta1);
    printf ("0-  X[2]: %g, Y[2]: %g\n",X[2],Y[2]);
    MPI_Bcast(X,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Y,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&theta0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&theta1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Recv(&ct, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf ("0-  sum: %g\n",ct);
    //ct=cost(N,X,Y,theta0,theta1);
    printf ("(%d) Theta; %g, %g  cost: %g\n",i,theta0,theta1,ct);
  }
  else{
    int i;
    double val;
    double sum = 0.0;
    double *vx = malloc(sizeof(double)* N);
    double *vy = malloc(sizeof(double)* N);
    double t0,t1;
    MPI_Bcast(vx,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(vy,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    printf ("1-  vx[2]: %g, vy[2]: %g\n",vx[2],vy[2]);
    MPI_Bcast(&t0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&t1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for(i=0;i<N;i++)
    {
      val = t0 + t1*vx[i] - vy[i];
      sum += val * val;
    }
    sum /= 2*N;
    MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    printf ("1-  sum: %g\n", sum);
  }  
  MPI_Finalize();
  return(0);
}
