#include "task_7.hpp"
using namespace std;

double f(double x) {
    //return sin( M_PI * x);
    //return x*exp(x);
    //return sin(M_PI * 1 * (x-1./2.) / (double)(10 - 1));
    //return ((exp(x)-exp(1.0))*sin(x));
    return (x*(1-x)*cos(x*x));
    //return x;
    //return x*(1 - x);
    //return (x*(1+x));
}

void printMatrix(const std::vector<double>& matrix) {
    int n = sqrt(matrix.size());
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << std::setw(10) << std::setprecision(4) << matrix[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(const std::vector<double>& vec){
    for (double val : vec) {
        std::cout << std::setw(10) << std::setprecision(4) << val << " ";
    }
    std::cout << std::endl;
}

void FillingNodes(std::vector<double>& xk){
    int N = xk.size(); //N+1
    N--;
    double h = 1./(double)(N-1);

    xk[0] = -h/2.;
    for (int i=1;i<N + 1;++i){
        xk[i] = xk[i-1] + h;
    }
}

void FillingBasicNodes(std::vector<double>& xk){
    int N = xk.size(); //N-1
    N++;
    double h = 1./(double)(N-1);
    xk[0] = h/2.;
    // X_1, .., X_N-1
    for (int i = 1; i < N - 1; ++i){
        xk[i] = xk[i - 1] + h;
    }
}

void FillingValues(std::vector<double>& xk, std::vector<double>& yk, std::function<double(double)> f){
    int N = xk.size(); //N+1
    N--;

    if (xk.empty()) {
        std::cerr<<"Please, fill nodes xk!"<<std::endl;
        throw std::runtime_error("xk not found");
    }
    for (int i=0;i<N+1;++i){
        yk[i] = f(xk[i]);
    }
}

void FillingBasicValues(std::vector<double>& xk, std::vector<double>& yk, std::function<double(double)> f){
    int N = xk.size(); //N-1
    //size_t N = xk.size();
    /*N++;

    if (xk.empty()) {
        std::cerr<<"Please, fill nodes xk!"<<std::endl;
        throw std::runtime_error("xk not found");
    }*/
    N++;
    for (int i=0;i<N-1;++i){
        yk[i] = f(xk[i]);
    }
}

// функция вычисляет базисную функцию
double Psi(int k, int n, int N){
    return sin(M_PI * n * (k-1./2.) / (double)(N - 1));
}

// функция вычисляет собственные значения матрицы А
double Lambda(int n, int N, double p){
    double lam = p - 2. * (double)(N-1) * (double)(N-1) * (cos(M_PI * n / (double)(N-1)) - 1.);
    return lam;
}

// функция вычисляет скалярное произведение fk на базисную функцию {Psi_k}^(n)
double ScalarProduct(int n, std::vector<double>& fk){
    int N = fk.size(); //N-1
    N++;
    double sp = 0.;
    double norm = 0.;

    for (int i = 0; i < N - 1; ++i){
        sp += fk[i] * Psi(i + 1, n, N) / ((double)(N-1));
        norm += Psi(i + 1, n, N) * Psi(i + 1, n, N) / ((double)(N-1));
    }
    sp = sp / (double)norm;
    return sp;
}

// функция реализует метод разложения вектора fk по базису Фурье для нахождения yk
// fk - правая часть уравнения
void Fourier(std::vector<double>& y, double p, std::vector<double>& fk){
    int N = fk.size(); //N-1
    N++;

    for (int k = 0; k < N - 1; ++k){
        y[k]=0.;
        for (int n = 1; n < N; ++n){
            y[k] += (ScalarProduct(n, fk) / Lambda(n, N, p)) * Psi(k + 1, n, N);
            //std::cout<<"with n = "<<n<<" coeff = "<<(ScalarProduct(n, fk) / Lambda(n, N, p) )<<std::endl;
        }
        //std::cout<<"with k = "<< k << ": y[k] - 1 = " << y[k] - 1. <<std::endl;
    }
    /*for (int k = 0; k < N - 1; ++k){
        std::cout<<"with k = "<< k << ": y[k] - 1 = " << y[k] - 1. <<std::endl;
    }*/
}

// функция реализует умножение матрицы на вектор
std::vector<double> MultiplyMatrixByVector(const std::vector<double>& matrix, const std::vector<double>& vec) {
    int n = sqrt(matrix.size());
    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += matrix[i*n+j] * vec[j];
        }
    }
    return result;
}

// функция заполняет матрицу А для задачи с фиксированным p
void BasisMatrixFill(double p, std::vector<double>& M){
    int N = sqrt(M.size());  //N-1
    N++;

    for (int i = 0; i < N-1; ++i) {
        for (int j = 0; j < N-1; ++j) {
            if (i == j){
                M[i * (N - 1) + j] = 2.0 * (double)((N-1) * (N-1)) + p + 500*sin(double(i));
                //M[i * (N - 1) + j] = 2.0 * (double)((N-1) * (N-1)) + p;
            }
            else if ((i - j == 1 || i - j == -1)){
                M[i * (N - 1) + j] = -(double)((N-1) * (N-1));
            }
            else{
                M[i * (N - 1) + j] = 0.0;
            }
        }
    }
    M[0] = 3.* (double)((N-1) * (N-1)) + p + 500*sin(0.);
    M[(N - 2) * (N - 1) + (N - 2)] = 3.* (double)((N-1) * (N-1)) + p + 500*sin(double(N - 2));
    /*M[0] = 3.* (double)((N-1) * (N-1)) + p;
    M[(N - 2) * (N - 1) + (N - 2)] = 3.* (double)((N-1) * (N-1)) + p;*/
}

// аналогично BasisMatrixFill, но у нас переменный коэффициент pk
void BasisMatrixFillWithVariableP(double p, std::vector<double>& M) {
    int N = sqrt(M.size()); //N-1
    N++;
    double mean = 0.;
    double pk;
    for (int i = 0; i < N-1; ++i) {
        pk = sin(M_PI * i / (double)(N - 1));
        //pk = 1e-5;
        //xpk = 50.;
        //pk = 1.;
        //pk = 0.;
        for (int j = 0; j < N-1; ++j) {
            if (i == j){
                M[i * (N - 1) + j] = 2.0 * (double)((N-1) * (N-1)) + 1.+pk*pk + p;
                //M[i * (N - 1) + j] = 2.0 * (double)((N-1) * (N-1)) + p + pk;
                //M[i * (N - 1) + j] = 2.0 * (double)((N-1) * (N-1)) + pk + 500*sin(double(i));;
            }
            else if ((i - j == 1 || i - j == -1)){
                M[i * (N - 1) + j] = -(double)((N-1) * (N-1));
            }
            else{
                M[i * (N - 1) + j] = 0.0;
            }
        }
        mean += 1+pk*pk;
        //mean += pk;
    }
    mean/=(double)(N-1);
    //pk = 1e-5;
    //pk = 50.;
    //pk = 1.;
    //pk = 0.;
    /*M[0] = 3.* (double)((N-1) * (N-1)) + p + pk;
    M[(N - 2) * (N - 1) + (N - 2)] = 3.* (double)((N-1) * (N-1)) + p + pk;*/
    /*M[0] = 3.* (double)((N-1) * (N-1)) + pk + 500*sin(0.);
    M[(N - 2) * (N - 1) + (N - 2)] = 3.* (double)((N-1) * (N-1)) + pk + 500*sin(double(N - 2));*/
    double pk_1 = 1. + sin(M_PI * 0 / (double)(N - 1));
    double pk_2 = 1. + sin(M_PI * (N- 2) / (double)(N - 1));
    M[0] = 3.* (double)((N-1) * (N-1)) + 1. + pk_1 * pk_1 + p;
    M[(N - 2) * (N - 1) + (N - 2)] = 3.* (double)((N-1) * (N-1)) + 1. + pk_2 * pk_2 + p;
    //return mean;
}

void FullMatrixFill(double p, std::vector<double>& M){
    int N = sqrt(M.size());  //N+1
    N--;

    for (int i = 1; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j){
                M[i * (N + 1) + j] = 2.0 * (double)((N-1) * (N-1)) + p;
                //M[i * (N + 1) + j] = 2.0 * (double)((N-1) * (N-1)) + p + 500*sin(double(i - 1));
            }
            else if ((i - j == 1 || i - j == -1)){
                M[i * (N + 1) + j] = -(double)((N-1) * (N-1));
            }
            else{
                M[i * (N + 1) + j] = 0.0;
            }
        }
    }
    // y_0 = - y_1
    M[0] = 1.;
    M[1] = 1.;
    M[1 * (N + 1)] = 0.;

    M[1 * (N + 1) + 1] = p + 3. * (double)((N-1) * (N-1));
    M[N * (N + 1) - 2] = p + 3. * (double)((N-1) * (N-1));

    /*M[1 * (N + 1) + 1] = 3.* (double)((N-1) * (N-1)) + p + 500*sin(double(1 - 1));
    M[(N - 1) * (N + 1) + N - 1] = 3.* (double)((N-1) * (N-1)) + p + 500*sin(double(N - 2));*/

    M[N * (N + 1) - 1] = 0.;

    //y_N = -y_{N - 1}
    M[N * (N + 1) + N - 1] = 1.;
    M[N * (N + 1) + N] = 1.;
}

void FullMatrixFillWithVariableP(std::vector<double>& M) {
    int N = sqrt(M.size()); //N+1
    N--;
    double pk;
    for (int i = 1; i < N; ++i) {
        //pk = sin(M_PI * i / (double)(N - 1));
        //pk = 1e-3;
        pk = 1.;
        //pk = 0.;
        for (int j = 0; j < N; ++j) {
            if (i == j){
                //M[i * (N + 1) + j] = 2.0 * (double)((N-1) * (N-1)) + 1.+pk*pk;
                M[i * (N + 1) + j] = 2.0 * (double)((N-1) * (N-1)) + pk;
                //M[i * (N + 1) + j] = 2.0 * (double)((N-1) * (N-1)) + pk + 500*sin(double(i - 1));
            }
            else if ((i - j == 1 || i - j == -1)){
                M[i * (N + 1) + j] = -(double)((N-1) * (N-1));
            }
            else{
                M[i * (N + 1) + j] = 0.0;
            }
        }
    }
    //pk = sin(M_PI * 0 / (double)(N - 1)) *  sin(M_PI * 0 / (double)(N - 1)) +1.;
    //pk = 1e-3;
    pk = 1.;
    //pk = 0.;

    // y_0 = - y_1
    M[0] = 1.;
    M[1] = 1.;
    M[1 * (N + 1)] = 0.;
    M[1 * (N + 1) + 1] = 3. * (double)((N-1) * (N-1)) + pk;
    M[N * (N + 1) - 2] = 3. * (double)((N-1) * (N-1)) + pk;

    /*M[1 * (N + 1) + 1] = 3. * (double)((N-1) * (N-1)) + pk + 500*sin(double(1 - 1));
    M[N * (N + 1) - 2] = 3. * (double)((N-1) * (N-1)) + pk + 500*sin(double(N - 1));*/

    M[N * (N + 1) - 1] = 0.;

    //y_N = -y_{N - 1}
    M[N * (N + 1) + N - 1] = 1.;
    M[N * (N + 1) + N] = 1.;
}

// функция BSolver реализует метод с предобуславливателем (использует Фурье для решения):
// умножает х на А, вычисляя r = b - Ах; разлагает r по базису psi^(n) (используя Фурье)
//  обновляет значения с помощью предобуславливателя
double BSolver( std::vector<double>& x, std::vector<double>& A,
                std::vector<double>& B, std::vector<double>& b,
                double tau, int mIter, std::vector<double>& mem,
                std::vector<double>& mem1, double p) {

    int N = x.size(); //N-1
    B[0]+=0;
    N++;

    for (int k = 0; k < N-1; ++k){
        x[k] = 0;
        mem[k] = 0;
        mem1[k] = 0;
    }
    for (int m = 0; m < mIter; ++m) {

        mem = MultiplyMatrixByVector(A, x);
        // b - Ax
        for (int j = 0; j < N-1; ++j) {
            mem[j] = b[j] - mem[j];
        }
        Fourier(mem1, p, mem);

        //  обновляет значения с помощью предобуславливателя
        for (int i = 0; i < N-1; ++i) {
            x[i] += tau * mem1[i];
        }
    }

    return ErNormInf(A, b, x, mem1);
}

// функция считает норму невязки r = b - Ax - это норма невязки в бесконечной норме (как максимум модуля)
double ErNormInf(std::vector<double>& A, std::vector<double>& b, std::vector<double>& x, std::vector<double>& mem){

    int N = x.size(); //N-1
    N++;
    double ans = 0.;
    mem = MultiplyMatrixByVector(A, x);

    for (int i = 0; i < N-1; i++){
        if(fabs((b[i] - mem[i])) > ans) ans = fabs((b[i] - mem[i]));
    }

    return ans;
}
