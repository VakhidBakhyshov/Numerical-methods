#include "task_2.hpp"

int pcalculate(int N){

double* xk      = new double[N+1];
double* fmemory = new double[N+1];
double* u       = new double[(N+1)*(N+1)];
double* d       = new double[(N+1)*(N+1)];
double* c       = new double[(N+1)*(N+1)];
double* phi     = new double[N+1];

double* xk2      = new double[2*N+1];
double* fmemory2 = new double[2*N+1];
double* u2       = new double[(2*N+1)*(2*N+1)];
double* d2       = new double[(2*N+1)*(2*N+1)];
double* c2       = new double[(2*N+1)*(2*N+1)];
double* phi2     = new double[2*N+1];

double h1 = 1/((double)N-1);
double h2 = 1/(double)((2*N)-1);

FillingNodes(xk, N);
FillingUMatrix(N, u, xk);
FillingDMatrix(N, d, u, phi);
FillingCMatrix(N, d, c, fmemory, phi);

//printMatrixCToConsole(c, N);

FillingNodes(xk2, 2*N);
FillingUMatrix(2*N, u2, xk2);
FillingDMatrix(2*N, d2, u2, phi2);
FillingCMatrix(2*N, d2, c2, fmemory2, phi2);

double err1 = normFunction(f, c, N);
std::cout << "err1 = " << fabs(err1) << std::endl;
double err2 = normFunction(f, c2, 2*N);
std::cout << "err2 = " << fabs(err2) << std::endl;

double a = log(err1/err2);
double b = log(h1/h2);
std::cout << "p =  " << fabs(a/b) << std::endl;

delete[] xk;
delete[] fmemory;
delete[] u;
delete[] d;
delete[] c;
delete[] phi;
delete[] xk2;
delete[] fmemory2;
delete[] u2;
delete[] d2;
delete[] c2;
delete[] phi2;

return 0;

}
