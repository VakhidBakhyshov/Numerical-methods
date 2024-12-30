#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <functional>

double f(double x);
void printMatrix(const std::vector<double>& matrix);
void printVector(const std::vector<double>& vec);

void FillingNodes(std::vector<double>& xk);
void FillingValues(std::vector<double>& xk, std::vector<double>& yk, std::function<double(double)> f);

void FillingBasicNodes(std::vector<double>& xk);
void FillingBasicValues(std::vector<double>& xk, std::vector<double>& yk, std::function<double(double)> f);

double Psi(int k, int n, int N);
double Lambda(int n, int N, double p);
double ScalarProduct(int n, std::vector<double>& fk);
void Fourier(std::vector<double>& y, double p, std::vector<double>& fk);


std::vector<double> MultiplyMatrixByVector(const std::vector<double>& matrix, const std::vector<double>& vec);

void BasisMatrixFill(double p, std::vector<double>& M);
void FullMatrixFill(double p, std::vector<double>& M);

void BasisMatrixFillWithVariableP(double p, std::vector<double>& M);
void FullMatrixFillWithVariableP(std::vector<double>& M);

double SearchQ(std::vector<double>& A);

double BSolver( std::vector<double>& x, std::vector<double>& A,
                std::vector<double>& B, std::vector<double>& b,
                double tau, int mIter, std::vector<double>& mem,
                std::vector<double>& mem1, double p);
double ErNormInf(std::vector<double>& A, std::vector<double>& b, std::vector<double>& x, std::vector<double>& mem);
