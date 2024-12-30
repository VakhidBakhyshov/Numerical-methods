#include "task_1.hpp"

double f(double x){
    //return ((exp(x)-exp(1.0))*sin(x));
    //return (x*(1-x)*cos(x*x));
    return (x*(1-x));
}


void Filling Nodes(double* xk, int N){
    double h = 1/(N-1.);
    for (int i = 0; i < N + 1; i++){
	xk[i] = -h/2. + i*h;
    }
}

void FillingValues(double* xk, double* yk, double (*f)(double), int N){
    for (int i=0;i<N+1;++i){
        yk[i] = f(xk[i]);
    }
}

double PhiCalculate(int n, int k, double h){
    return sin(M_PI * n * (-h/2. + k*h));
}

void PhiVectorCalculate(int N, int n, double* phi){  //phi_k^{(n)} for fixed k
    double h = 1/(N-1.);
    for (int k=1;k<=N;++k){
        phi[k] = PhiCalculate(n, k, h);
    }
}

double Scalar(double* ar1,double* ar2, int N){
    double s1 = 0., s2 = 0.;
    for (int k=1; k<N; ++k){
        s1 += ar1[k] * ar2[k];
        s2 += ar1[k] * ar1[k];
    }
    s1 = s1/double(s2);
    return s1;
}

void CoeffCalculate(int N, double* yk, double* phi, double* cn){
    for (int n=1; n<=N; ++n){
        PhiVectorCalculate(N, n, phi);
        cn[n] = Scalar(phi, yk, N);
    }
}

double FourierCompute(double* cn, int N, double x){
    double res = 0;
    for(int n = 1; n <= N; ++n){
        res += cn[n] * sin( M_PI * n * x);
    }
    return res;
}

void WriteToFile(const std::string& filename, int N, double* xk, double* yk, double* cn, double* phi) {
    double h = 1/(N-1.);
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        outFile << "         xk                   yk                      yk*              "<<std::endl;
        for (int i = 1; i < N; ++i){
        //for (int i = 1; i < N + 1; ++i){
            CoeffCalculate(N, yk, phi, cn);
            FourierCompute(cn, N, (- h/2. + i* h));
            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xk[i] << " " \

            << std::setw(20) << yk[i] << " " \
            << std::setw(20) << FourierCompute(cn, N, xk[i]) << std::endl;
        }
        outFile.close();
        std::cout << "Information successfully written to file!" << std::endl;
    }
    else {
        std::cerr << "Error opening file" << std::endl;
    }
}

double NormFunction(double (*f)(double), double* cn, int N){
    //double h = 1/100.;
    double h = 1/(N - 1.);
    double max = -1;
    double delta = 0;
    for (double x = 0.; x < 1.; x += h){
        delta = fabs(f(x) - FourierCompute(cn,N,x));
        if (delta > max) max = delta;
    }
    return max;
}
