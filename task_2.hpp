#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

double AnalyticalSolution(double coefficient, double init_value, double x_coord);
double ExplicitEulerMethod(double initial_value, double coefficient, int steps);
double ImplicitEulerMethod(double initial_condition, double param, int divisions);
double TrapezoidalApproach(double y_init, double A_param, int N_intervals);
double LeapfrogTechnique(double y0, double alpha, int num_steps);