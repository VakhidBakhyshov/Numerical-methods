#include "task_3.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <fstream>


void Filling_x_Grid(std::vector<double>& x_k, double h){ //x_Grid
    int M = x_k.size();
    for (int j = 0; j < M; ++j) {
        x_k[j] = -h/2. + j*h;
    }
}

void Filling_t_Grid(std::vector<double>& t_k, double tau){ //t_Grid
    int N = t_k.size();
    for (int i = 0; i < N; ++i) {
        t_k[i] = i * tau;
    }
}

double func(double t, double x) {//f(t, x)

    //return exp(-M_PI * M_PI * t) * (- M_PI * M_PI * x * (x - 1.0) - 2); //func for 1st case
    //return 0.0; //func for 2nd case
    return exp(-M_PI * M_PI * t) * ((1 - M_PI * M_PI) * (x - 1) * sin(x) - 2 * cos(x)); //func for 3rd case
}

double initialValue(double x){ //u0 = u(x, 0)

    //return x * (x - 1.0); //u0 for 1st case
    //return sin(M_PI * x); //u0 for 2nd case
    return sin(x) * (x - 1); //u0 for 3rd case
}

double analyticalSolution(double t, double x){ //u_{exact}

    //return exp(-M_PI * M_PI * t) * x * (x - 1.0); //u_{exact} for 1st case
    //return exp(-t * M_PI * M_PI) * sin(M_PI * x); //u_{exact} for 2nd case
    return exp(-t *M_PI * M_PI) * sin(x) * (x - 1); //u_{exact} for 3rd case
}

void FillingInitialZeroLayer(int M, double h, std::vector<double>& u){
    for (int j = 0; j <= M; ++j) {
        u[0*(M+1)+j] = initialValue(-h/2. + j*h); //u0(-h/2. + j*h)
    }
}

void ExplicitScheme(int N, int M, double tau, double h, std::vector<double>& u){
    FillingInitialZeroLayer(M, h, u);

    double ratio = tau / (h * h);

    // The explicit scheme is stable when tay <= h^2 / 2, otherwise when tau > h^2 / 2 - explicit scheme is not stable, i.e. unstable
    if (ratio > 0.5){
        std::cerr << "Explicit scheme is not stable, because tau / h^2 = " << ratio << " > 0.5)\n";
    }

    for (int i = 0; i < N; ++i){
        for (int j = 1; j < M; ++j){
            u[(i+1)*(M+1)+j] = u[i*(M+1)+j] + ratio * (u[i*(M+1)+j-1] - 2.0 * u[i*(M+1)+j] + u[i*(M+1)+j+1]) + tau * func(i * tau, -h/2. + j*h);
        }

        // 2 boundary conditions: u_0 + u_1 = 0, u_{M-1} + u_M = 0
        u[(i+1)*(M+1)+0] = - u[(i+1)*(M+1)+1];
        u[(i+1)*(M+1)+M] = - u[(i+1)*(M+1)+M-1];

        /*if(i < N){
            for (int j=0;j<M+1;j++) std::cout <<u[(i+1)*(M+1)+j] <<" ";
            std::cout << std::endl;
            for (int j=0;j<M+1;j++) std::cout << analyticalSolution((i+1)*tau,-h/2.0+j*h) << " ";
            std::cout << std::endl << std::endl;
        }*/
    }
}


void setupTridiagonalSystem(std::vector<double>& upper, std::vector<double>& middle, std::vector<double>& lower, int M, double h, double tau){

    double coeff_1 =  -1.0 / (h * h);
    double coeff_2 = 1.0 / tau + 2.0 / (h * h);

    for (int i = 1; i < M; ++i){
        lower[i] = coeff_1;
        middle[i]  = coeff_2;
        upper[i] = coeff_1;
    }

    //1st boundary condition: u_0 + u_1 = 0
    lower[0] = 0.0;
    middle[0] = 1.0;
    upper[0]  = 1.0;

    //2nd boundary condition: u_{M-1} + u_M = 0
    lower[M] = 1.0;
    middle[M] = 1.0;
    upper[M]= 0.0;
}


void MethodProgonki(std::vector<double>& upper, std::vector<double>& middle, std::vector<double>& lower, std::vector<double>& y, int M){

    std::vector<double> a(M+1);
    std::vector<double> b(M+1);

    for (int i = 0; i < M; ++i){
        if (i == 0){
                a[0] = -upper[0] / middle[0];
                b[0] = y[0] / middle[0];
        }
        else{
                double coeff = middle[i] + lower[i] * a[i - 1];
                a[i] = -upper[i] / coeff;
                b[i] = (y[i] - lower[i] * b[i - 1]) / coeff;
        }
    }

    for (int i = M; i >= 0; --i){
        if (i == M) y[M] = (y[M] - lower[M] * b[M - 1]) / (middle[M] + lower[M] * a[M - 1]);
        else y[i] = a[i] * y[i + 1] + b[i];
    }
}


void ImplicitScheme(int N, int M, double tau, double h, std::vector<double>& u){

    FillingInitialZeroLayer(M, h, u);

    std::vector<double> f(M + 1);
    std::vector<double> upper(M + 1);
    std::vector<double> middle(M + 1);
    std::vector<double> lower(M + 1);

    setupTridiagonalSystem(upper, middle, lower, M, h, tau);

    //double coeff = (h * h) / 8;

    for (int i = 1; i <= N; ++i){
        //1st boundary condition: u_0 + u_1 = 0 -> f_0 = 0
        f[0] = 0.0;
        //2nd boundary condition: u_{M-1} + u_M = 0 -> f_M = 0
        f[M] = 0.0;

        //f[0] = coeff * ((u[i*(M+1) + 0] - u[(i-1)*(M+1) + 0]) / tau - func((i-1) * tau, 0));
        //f[M] = coeff * ((u[i*(M+1) + M] - u[(i-1)*(M+1) + M]) / tau - func((i-1) * tau, 1));

        for (int j = 1; j < M; ++j){
            f[j] = u[(i-1)*(M+1)+j] / tau + func(i * tau, -h/2.0 + j*h);
        }

        MethodProgonki(upper, middle, lower, f, M);

        for (int j = 0; j <= M; ++j) u[i*(M+1)+j] = f[j];
    }
}


double norma(std::vector<double>& u, int N, int M, double h, double T){
    double norma = 0.0;

    for (int i = 0; i <= M; ++i){
        double analytical = analyticalSolution(T,-h/2. + i*h);
        double diff = u[N*(M+1)+i] - analytical; //on N layer we calculate error
        norma += diff * diff;
    }

    return sqrt(norma * h);
}


void convergence_rate(double T, const int type_of_scheme, int key){
    std::ofstream file_4("rate1.txt");
    if (!file_4.is_open()) throw std::runtime_error("Error in opening rate.txt\n");

    std::ofstream file_5("rate2.txt");
    if (!file_5.is_open()) throw std::runtime_error("Error in opening rate.txt\n");

    if (type_of_scheme == 1){
        file_4 << std::fixed
            << std::setw(8) << "N" << std::setw(8) << "M" << std::setw(15) << "tau"
            << std::setw(15) << "h" << std::setw(15) << "Norma" << std::endl;
    }
    else if (type_of_scheme == 2){
        file_5 << std::fixed
            << std::setw(8) << "N" << std::setw(8) << "M" << std::setw(15) << "tau"
            << std::setw(15) << "h" << std::setw(15) << "Norma" << std::endl;
    }

    int number = 6;
    std::vector<double> array_M(number);
    for (int i = 0; i < number; i++) array_M[i] = 5 * pow(2, i);
    std::vector<double> array_N(number);

    if (key == 1){ //for explicit scheme we have M -> 2M, N -> 2N, i.e. M,N -> 2M,2N
        for (int i = 0; i < number; ++i) array_N[i] = 50000 * pow(2, i);
    }
    else if (key == 2){ //for implicit scheme we have M -> 2M, N -> 2N, i.e. M, N -> 2M,4N
        for (int i = 0; i < number; ++i) array_N[i] = 2 * pow(array_M[i], 2);
    }
    else if (key == 3){ //for impiciy scheme we have M = N
        for (int i = 0; i < number; ++i){
                //array_N[i] = 50000 * pow(2, i);
                array_N[i] = array_M[i];
        }
    }
    else std::cout << "error in choosing right case for calculating convergence rate\n" << std::endl;

    for (int i = 0; i < number; ++i){

        int N = array_N[i];
        int M = array_M[i];

        double tau = T / N;
        double h = 1.0 / (double)(M - 1);

        std::vector<double> u((N+1)*(M+1));

        if (type_of_scheme == 1) ExplicitScheme(N, M, tau, h, u);
        else if (type_of_scheme == 2) ImplicitScheme(N, M, tau, h, u);
        else std::cout << "unexisted type of scheme. enter the right type of scheme\n" << std::endl;

        double l2_norma = norma(u, N, M, h, T);

        if (type_of_scheme == 1){
                file_4 << std::scientific << std::setfill(' ')
                        << std::setw(8) << N << std::setw(8) << M << std::setw(15) << tau
                        << std::setw(15) << h << std::setw(15) << l2_norma << std::endl;
        }
        else if (type_of_scheme == 2){
                file_5 << std::scientific << std::setfill(' ')
                        << std::setw(8) << N << std::setw(8) << M << std::setw(15) << tau
                        << std::setw(15) << h << std::setw(15) << l2_norma << std::endl;
        }

    }

    if (type_of_scheme == 1) std::cout << "Statistics about convergence rate succesfully written to rate1.txt\n";
    else if (type_of_scheme == 2) std::cout << "Statistics about convergence rate succesfully written to rate2.txt\n";
    file_4.close();
    file_5.close();
}