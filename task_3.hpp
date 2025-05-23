#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>

void Filling_x_Grid(std::vector<double>& x_k, double h);
void Filling_t_Grid(std::vector<double>& t_k, double tau);

double func(double t, double x);
double initialValue(double x);
double analyticalSolution(double t, double x);
void FillingInitialZeroLayer(int M, double h, std::vector<double>& u);

void ExplicitScheme(int N, int M, double tau, double h, std::vector<double>& u);

void setupTridiagonalSystem(std::vector<double>& upper, std::vector<double>& middle, std::vector<double>& lower, int M, double h, double tau);
void MethodProgonki(std::vector<double>& upper, std::vector<double>& middle, std::vector<double>& lower, std::vector<double>& y, int M);

void ImplicitScheme(int N, int M, double tau, double h, std::vector<double>& u);

double norma(std::vector<double>& u, int N, int M, double h, double T);
void convergence_rate(double T, const int type_of_scheme, int key);