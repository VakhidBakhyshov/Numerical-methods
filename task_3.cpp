#include "task_3.hpp"

double f(double x){
    //return((exp(x)-exp(1.0))*sin(x));
    //return x;
    //return (x*x+sin(x))*cos(3*x);
    //return 1./(1.+25.*x*x); //Runge's function
    return std::fabs(x);
    //return pow (x , 10) + 5 * pow( x , 8) -2 * pow(x , 6) + 3 * pow (x , 5) + 2 * pow (x , 3) + x * x + 11;
}

int GenerateEquidistantNodes(double a, double b, double (*f)(double), std::vector<double>& nodes, std::vector<double>& values) {
    int N = nodes.size();
    /*if (N < 2 || a >= b) {
        return -1;
    }*/
    double step = (b - a) / double(N - 1);
    double start = a;

    for (int i = 0; i < N; ++i) {
        nodes[i] = start;
        values[i] = f(start);
        start += step;
    }
    return 0;
}

int GenerateChebyshevNodes(double a, double b, double (*f)(double), std::vector<double>& nodes, std::vector<double>& values){
    int N = nodes.size();
    /*if (N < 1 || a >= b) {
        return -1;
    }*/

    double a1 = (b + a) / 2, a2 = (b - a) / 2; //[a,b] <-> [-1,1]
    for (int i = 0; i < N; ++i) {
        nodes[N - 1 - i] = a1 + a2 * cos((2 * i + 1) * M_PI / (2 * N));
        values[N - 1 - i] = f(nodes[N - 1 - i]);
    }
    return 0;
}

int MatrixFill(std::vector<double>& matrix, const std::vector<double>& nodes){
    int n = sqrt(matrix.size());
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i*n+j] = pow(nodes[i],j);
        }
    }
    return 0;
}

int VecFill(std::vector<double>& vec, const std::vector<double>& values){
    int n = vec.size();
    for (int i = 0; i < n; ++i) {
        vec[i] = values[i];
    }
    return 0;
}

double PnCalculation(std::vector<double>& coeffs, double x){
    int n = coeffs.size();
    double ans = 0.;
     for (int i = 0; i < n; ++i) {
            ans += coeffs[i] * pow(x,i);
    }
    return ans;
}

double LnCalculation(const std::vector<double>& x, const std::vector<double>& y, double x_value) {
    int n = x.size();
    double ans = 0.0;
    for (int i = 0; i < n; ++i) {
        double tmp = y[i];
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                tmp *= (x_value - x[j]) / (x[i] - x[j]);
            }
        }
        ans += tmp;
    }
    return ans;
}


void  WriteToFilePn(double a, double b, const std::string& filename, std::vector<double>& coeffs){
    std::ofstream outFile(filename);
    //if (outFile.is_open()) {
        int N = coeffs.size();
        double h = (b - a) / double(N - 1);
        h /= 3.;
        double xi = a-h, tmp1,tmp2,err;

        outFile<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
        <<std::setw(9)<<" "<<"f(x)"<<std::setw(9)<<" "\
        <<std::setw(9)<<" "<<"Pn"<<std::setw(9)<<" "\
        <<std::setw(9)<<" "<<"err"<<std::setw(9)<<" "<<std::endl;

        int M = (b-a)/(h*3);
        for (int i = 0; i < M; ++i){
            xi += h;
            tmp1 = f(xi);
            tmp2 =  PnCalculation(coeffs,xi);
            err = std::fabs(tmp1-tmp2);

            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;

            xi += h;
            tmp1 = f(xi);
            tmp2 =  PnCalculation(coeffs,xi);
            err = std::fabs(tmp1-tmp2);

            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;

            xi += h;
            tmp1 = f(xi);
            tmp2 =  PnCalculation(coeffs,xi);
            err = std::fabs(tmp1-tmp2);

            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;

        }
    //}
    //else {
    //    std::cerr << "Error opening file" << std::endl;
    //}
}



void  WriteToFileLn(double a, double b, const std::string& filename, std::vector<double>& nodes, std::vector<double>& values){
    std::ofstream outFile(filename);
    //if (outFile.is_open()) {
        int N = values.size();
        double h = (b - a) / double(N - 1);
        h /= 3.;
        double xi = a-h, tmp1,tmp2,err;

        outFile<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
        <<std::setw(9)<<" "<<"f(x)"<<std::setw(9)<<" "\
        <<std::setw(9)<<" "<<"Ln"<<std::setw(9)<<" "\
        <<std::setw(9)<<" "<<"err"<<std::setw(9)<<" "<<std::endl;

        int M = (b-a)/(h*3);
        for (int i = 0; i < M; ++i){
            xi += h;
            tmp1 = f(xi);
            tmp2 = LnCalculation(nodes, values, xi);
            err = std::fabs(tmp1-tmp2);

            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;

            xi += h;
            tmp1 = f(xi);
            tmp2 =  LnCalculation(nodes, values, xi);
            err = std::fabs(tmp1-tmp2);

            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;

            xi += h;
            tmp1 = f(xi);
            tmp2 =  LnCalculation(nodes, values, xi);
            err = std::fabs(tmp1-tmp2);

            outFile << std::setprecision(15) << std::fixed \
            << std::setw(20) << xi   << " " \
            << std::setw(20) << tmp1 << " " \
            << std::setw(20) << tmp2 << " " \
            << std::setw(20) << err  << std::endl;

        }
    //}
    //else {
    //    std::cerr << "Error opening file" << std::endl;
    //}
}
