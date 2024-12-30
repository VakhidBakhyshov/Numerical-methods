#ifndef SOLUTION_H
#define SOLUTION_H
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>


# define M_PI 3.14159265358979323846

//double f(double x, double y, int k);
//double trueValue(double xa, double xb, double ya, double yb, int k);

void triangulationFixed(double xa, double xb, double ya, double yb, int N);
double f(double x, double y, int p, int q);
double trueValue(double xa, double xb, double ya, double yb, int p, int q);
void triangulation(double xa, double xb, double ya, double yb, int N);
void readingQuadratureFile(double *CoordsAndWeights);
void getVertexCoordinates(double xa, double xb, double ya, double yb, int m, int N, double &xm, double &ym);

//double integrateAccordingToTheQuadratureOfTheCenterOfMass(double (*f)(double, double, int), double xa, double xb, double ya, double yb, int k, int N, double *VertexCoords);
double integrateAccordingToTheQuadratureOfTheCenterOfMass(double (*f)(double, double, int, int),
                                                            double xa, double xb, double ya, double yb,
                                                            int N, double *VertexCoords, int p, int q);

void searchAffineTransformation(double* VertexCoords, double *CoefOfTransform);
//double integrateAccordingToTheQuadratureFromTheFile(double (*f)(double, double, int), double xa, double xb, double ya, double yb, int k, int N, double *VertexCoords, double *CoefOfTransform, double *CoordsAndWeights);
double integrateAccordingToTheQuadratureFromTheFile(double (*f)(double, double, int, int), double xa, double xb, double ya, double yb, int N, double *VertexCoords, double *CoefOfTransform, double *CoordsAndWeights, int p, int q);

//void writeToFileForP(double (*f)(double, double, int), double xa, double xb, double ya, double yb, int k, int N, int numTests, double *VertexCoords, double *CoefOfTransform, double *CoordsAndWeights, int typeOfQuadrature);
void writeToFileForP(double (*f)(double, double, int, int),
                        double xa, double xb, double ya, double yb, int N,
                        int numTests, double *VertexCoords, double *CoefOfTransform,
                        double *CoordsAndWeights, int typeOfQuadrature, int p, int q);

#endif
