#include "task_2.hpp"

double f(double x, double y){
    //return x * (1 - x) * y * (1 - y);
    //return x * (x - 1) * y * (1 - y);
    //return x * (1 - x) * y * (1 - y) * cos(x * x) * cos(y * y);
    return sin(5 * M_PI * x) * sin(M_PI * y);
    //return sin(2 * M_PI * x) * sin(3 * M_PI * y);
}

void FillingNodes(double* xk, int N){
    double h = 1/((double)N-1.);
    for (int i = 0; i < N+1; i++){
        xk[i] = -h/2 + i*h;
    }
}

void FillingUMatrix(int N, double* U, double* xk){
    for (int i = 0; i < N+1; i++){
        for (int j = 0; j < N+1; j++){
            U[i * (N+1) + j] = f(xk[i], xk[j]);
        }
    }
}

double PhiCalculate(int n, int k, double h){
    return sin(M_PI * n * (-h/2. + k*h));
}

void PhiVectorCalculate(int N, int n, double* phi){  //phi_k^{(n)} for fixed n
    double h = 1/((double)N-1.);
    for (int k = 0; k < N + 1; k++){
        phi[k] = PhiCalculate(n, k, h);
    }
}

double Scalar(double* ar1,double* ar2, int N){
    double h = 1/((double)N-1.), s1 = 0., s2 = 0.;
    for (int i = 1;i < N; i++){
        s1 += ar1[i] * ar2[i] * h;
        s2 += ar1[i] * ar1[i] * h;
    }
    s1 = s1 / (double)s2;
    return s1;
}

void CoeffCalculate(int N, double* yk, double* phi, double* cn){
    for (int n = 1; n < N; n++){
        PhiVectorCalculate(N, n, phi);
        cn[n] = Scalar(phi, yk, N);
    }
}

void FillingDMatrix(int N, double* D, double* U, double* phi){
    for (int i = 0; i < N + 1; i++){
        D[i] = 0;
        D[(N + 1) * N + i] = 0;
    }
    for (int i = 1; i < N; i++){
        CoeffCalculate(N, U + i * (N+1), phi, D + i * (N+1));
    }
}

void FillingCMatrix(int N, double* D, double* C, double* fmemory, double* phi){
    for (int i = 0; i < N + 1; i++){
        C[i] = 0;
        C[(N + 1) * N + i] = 0;
    }
    for (int i = 1; i < N; i++){
        for (int j = 0; j < N+1; j++){
            fmemory[j] = D[j*(N+1) + i];
        }
        CoeffCalculate(N, fmemory, phi, C + i * (N+1));
    }
}

double Calc2DFourier(double* C, int N, double x, double y){
    double res = 0;
    for (int i = 1; i < N; i++){
        for (int j = 1; j < N; j++){
            res += C[i * (N+1) + j] * sin(M_PI * j * x) * sin(M_PI * i * y);
        }
    }
    return res;
}

double normFunction(double (*f)(double, double), double *C, int N){
    //double h = 1./100.;
    double h = 1/((double)N-1.);
    double max = -1.;
    double delta = 0.;
    for (double x = 0.; x < 1.; x += h){
        for (double y = 0.; y < 1.; y += h){
            delta = fabs(f(x, y) - Calc2DFourier(C, N, x, y));
            if (delta > max){
                max = delta;
            }
            //std::cout << "x = " << x  << ", " << "y = " << y << ", " << "delta = " << delta << "\n" << std::endl;
            //std::cout << "max = " << max << "\n" << std::endl;
        }
    }
    return max;
}

void  WriteToFile(const std::string& filename, int N, double* xk, double* U, double* C, double* D, double* fmemory, double* phi){
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        //double h = 1/(N-1.);
        double xi = 0, yi = 0;
        double deltax = 0, deltay = 0;
        outFile<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
        <<std::setw(10)<<" "<<"y"<<std::setw(10)<<" "\
        <<std::setw(9)<<" "<<"f(x,y)"<<std::setw(9)<<" "\
        <<std::setw(6)<<" "<<"Fourier "<<std::setw(6)<<" "<<std::endl;
        FillingNodes(xk, N);
        FillingUMatrix(N, U, xk);
        FillingDMatrix(N, D, U,phi);
        FillingCMatrix(N, D, C, fmemory, phi);
        for (int i = 1; i < N-1; ++i){
            for (int j = 1; j < N; ++j){
                xi = xk[i];
                deltax = xk[i+1] - xi;
                deltax /= 2.;
                yi = xk[j];
                deltay = xk[j+1] - yi;
                deltay /= 2.;
                outFile << std::setprecision(15) << std::fixed \
                << std::setw(20) << xi << " " \
                << std::setw(20) << yi << " " \
                << std::setw(20) << f(xi,yi) << " " \
                << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

                xi += deltax;
                yi += deltay;
                outFile<< std::setprecision(15) << std::fixed \
                << std::setw(20) << xi << " " \
                << std::setw(20) << yi << " " \
                << std::setw(20) << f(xi,yi) << " " \
                << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

                xi += deltax;
                yi += deltay;
                outFile<< std::setprecision(15) << std::fixed \
                << std::setw(20) << xi << " " \
                << std::setw(20) << yi << " " \
                << std::setw(20) << f(xi,yi) << " " \
                << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;
            }
        }
    }
    else {
        std::cerr << "Error opening file" << std::endl;
    }
}

double WriteToConsole(int N, double* xk, double* U, double* C, double* D, double* fmemory, double* phi){
    //double h = 1/((double)N-1.);
    //double xi = 0, yi = 0;
    //double deltax = 0, deltay = 0;
    std::cout<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
    <<std::setw(10)<<" "<<"y"<<std::setw(10)<<" "\
    <<std::setw(9)<<" "<<"f(x,y)"<<std::setw(9)<<" "\
    <<std::setw(6)<<" "<<"Fourier "<<std::setw(6)<<" "<<std::endl;

    clock_t start=clock();
    FillingNodes(xk, N);
    FillingUMatrix(N, U, xk);
    FillingDMatrix(N, D, U,phi);
    FillingCMatrix(N, D, C, fmemory, phi);
    clock_t end=clock();
    double duration =(double)(end-start)/CLOCKS_PER_SEC;

    /*for (int i = 1; i < N-1; ++i){
        for (int j = 1; j < N; ++j){
            xi = xk[i];
            deltax = xk[i+1] - xi;
            deltax /= 2.;
            yi = xk[j];
            deltay = xk[j+1] - yi;
            deltay /= 2.;
            std::cout << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << yi << " " \
            << std::setw(20) << f(xi,yi) << " " \
            << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

            xi += deltax;
            yi += deltay;
            std::cout << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << yi << " " \
            << std::setw(20) << f(xi,yi) << " " \
            << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;

            xi += deltax;
            yi += deltay;
            std::cout << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi << " " \
            << std::setw(20) << yi << " " \
            << std::setw(20) << f(xi,yi) << " " \
            << std::setw(20) << Calc2DFourier(C, N, xi, yi) << std::endl;
        }
    }*/
    std::cout<<std::endl;
    return duration;
}

void printMatrixCToConsole(double *C, int N){
    //std::cout << "========================================================================= " << " Matrix C " << " =========================================================================" << std::endl;
    for (int i = 0; i < N+1; i++){
        for (int j = 0; j < N+1; j++){
            if (i == 0 || j == 0 || i == N || j == N){
                std::cout << std::setw(12)  << 0 << " ";
            }
            else{
                std::cout << std::setw(12) << std::setprecision(3) << C[i * (N + 1) + j] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

