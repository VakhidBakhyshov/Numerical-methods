//#ifndef SOLUTION_H
//#define SOLUTION_H
#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string>


# define M_PI 3.14159265358979323846
#define MAX_ITERATIONS 1500

double f(double x);
void fillingNodes(double a, double b, std::vector<double>& Nodes, int type);
void fillingValuesInNodes(std::vector<double>& Nodes, std::vector<double>& ValuesInNodes, int N);
void bubbleSort(std::vector<double>& Nodes, int N);
void fillingMatrixA(std::vector<double>& A, std::vector<double>& Nodes, int N);
void printMatrixAToConsole(std::vector<double>& A, int N);
void fillingVectorB(std::vector<double>& B, std::vector<double>& ValuesInNodes, int N);
void printCoefToConsole(std::vector<double>& Coef, int N);
double pnCalculating(std::vector<double>& Coef, int N, double x);
void writeToFilePn(std::vector<double>& Coef, int N, std::vector<double>& ShreddedGrid, int altcount);
void fillingSigma(std::vector<double>& sigma, std::vector<double>& Nodes, int altcount, int N);
void fillingValuesInSigma(std::vector<double>& sigma, std::vector<double>& ValuesInSigma, int altcount);
bool maxDeviation(std::vector<double>& Nodes, std::vector<double>& sigma, std::vector<double>& Coef, std::vector<double>& ValuesInNodes, std::vector<double>& ValuesInSigma, int N, int altcount);
void fillingShreddedGrid(std::vector<double>& ShreddedGrid, std::vector<double>& Nodes, int N);
void fillingValuesInShreddedGrid(std::vector<double>& ValuesInShreddedGrid, std::vector<double>& ShreddedGrid, int N);
void printVectorToConsole(std::vector<double>& Vector, int N);
//double CalcPolynom(std::vector<double>& Coef, double x, int N);

//#endif
