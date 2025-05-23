#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <functional>

double f(double x);

void FillingNodes(std::vector<double>& xk);
void ComputeGridPoints(std::vector<double>& grid_nodes);

void FillingValues(std::vector<double>& xk, std::vector<double>& yk, std::function<double(double)> f);
void CalculateFunctionValues(const std::vector<double>& grid_nodes,
                           std::vector<double>& function_values,
                           const std::function<double(double)>& math_func);

double Psi(int k, int n, int N);
double Eigenvalues(int n, int N, double p);

double ScalarProduct(int n, const std::vector<double>& func_values);

void MethodFourier(std::vector<double>& result, double p,
                          const std::vector<double>& func_values);

void InitializeMatrix(std::vector<double>& matrix, double p);

std::vector<double> MatrixVectorProduct(const std::vector<double>& mat,
                                      const std::vector<double>& vec);

double ComputeResidualNorm(const std::vector<double>& matrix,
                         const std::vector<double>& right_hand_side,
                         const std::vector<double>& solution_vector_x,
                         int size,
                         std::vector<double>& tmp_vector);

double FillTridiagonalMatrixWithVariableParameter(std::vector<double>& matrix);