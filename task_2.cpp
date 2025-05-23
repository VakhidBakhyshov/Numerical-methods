#include "task_2.hpp"

double AnalyticalSolution(double coefficient, double init_value, double x_coord) {
    return init_value * exp(-coefficient * init_value * x_coord);
}

double ExplicitEulerMethod(double initial_value, double coefficient, int steps) {
    double step_size = 1.0 / steps;
    double current_value = initial_value;
    double error_accumulator = 0.0;

    for (int iteration = 1; iteration < steps; ++iteration) {
        current_value *= (1 - coefficient * step_size);

        double exact_solution = AnalyticalSolution(coefficient, initial_value,
                                                        iteration * step_size);
        error_accumulator += pow(exact_solution - current_value, 2);
    }

    return sqrt(error_accumulator * step_size);
}

double ImplicitEulerMethod(double initial_condition, double param, int divisions) {
    double delta = 1.0 / divisions;
    double solution = initial_condition;
    double squared_error = 0.0;

    for (int step = 1; step < divisions; ++step) {
        solution /= (1.0 + param * delta);

        double theoretical_value = AnalyticalSolution(param, initial_condition,
                                                          step * delta);
        squared_error += pow(theoretical_value - solution, 2);
    }

    return sqrt(squared_error * delta);
}

double TrapezoidalApproach(double y_init, double A_param, int N_intervals) {
    double h_step = 1.0 / N_intervals;
    double numerical_sol = y_init;
    double l2_norm = 0.0;

    for (int idx = 1; idx < N_intervals; ++idx) {
        numerical_sol *= (2.0 - A_param * h_step) / (2.0 + A_param * h_step);

        double exact_val = AnalyticalSolution(A_param, y_init, idx * h_step);
        l2_norm += pow(exact_val - numerical_sol, 2);
    }

    return sqrt(l2_norm * h_step);
}

double LeapfrogTechnique(double y0, double alpha, int num_steps) {
    double dt = 1 / num_steps;
    double prev_solution = y0;
    double curr_solution = y0 * (1 - alpha * dt);
    double error_metric = 0.0;

    for (int k = 2; k < num_steps; ++k) {
        double next_solution = prev_solution - 2.0 * alpha * dt * curr_solution;
        prev_solution = curr_solution;
        curr_solution = next_solution;

        double true_solution = AnalyticalSolution(alpha, y0, k * dt);
        error_metric += (true_solution - curr_solution) * (true_solution - curr_solution);
    }

    return sqrt(error_metric * dt);
}
