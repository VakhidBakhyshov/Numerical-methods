#include "task_1.hpp"

double f(double x) {
    //return sin( M_PI * x);
    //return x*x*(1-x)*exp(x);
    //return x+1-x;
    return 1 - ((exp(x) + exp(1-x))/(1+exp(1)));
    //return sin(M_PI * 1 * (x-1./2.) / (double)(10-1));
    //return ((exp(x)-exp(1.0))*sin(x));
    //return (x*(1-x)*cos(x*x));
}

void ComputeGridPoints(std::vector<double>& grid_nodes){
    int num_intervals = grid_nodes.size();
    double h = 1.0 / double(num_intervals);

    grid_nodes[0] = - h/2.0;

    for (int k = 1; k < num_intervals; ++k) {
        grid_nodes[k] = grid_nodes[k - 1] + h;
    }
}

void CalculateFunctionValues(const std::vector<double>& grid_nodes,
                           std::vector<double>& function_values,
                           const std::function<double(double)>& f) {

    int num_points = grid_nodes.size();
    for (int k = 0; k < num_points; ++k) {
        function_values[k] = f(grid_nodes[k]);
    }
}

double Psi(int k, int n, int N){
    return sin(M_PI * n * (k - 0.5) / (double)(N-1));
}

double Eigenvalues(int n, int N, double p){
    return p - 2. * (double)(N-1) * (double)(N-1) * (cos(M_PI * n / (double)(N-1)) - 1.); // lambda, i.e. eigenvalue
}

double ScalarProduct(int n, const std::vector<double>& func_values) {
    int num_points = func_values.size() + 1;
    double scalar_product = 0.0;

    for (int i = 0; i < num_points - 1; ++i) {
        scalar_product += 2. * func_values[i] *
                     Psi(i + 1, n, num_points) /
                     (double)(num_points - 1);
    }
    return scalar_product;
}

void MethodFourier(std::vector<double>& result, double p,
                          const std::vector<double>& func_values) {

    int num_terms = func_values.size() + 1;
    result.resize(num_terms - 1);

    for (int k = 0; k < num_terms - 1; ++k) {
        result[k] = 0.0;
        for (int n = 1; n < num_terms; ++n) {
            double coeff = ScalarProduct(n, func_values) /
                          Eigenvalues(n, num_terms, p);
            result[k] += coeff * Psi(k + 1, n, num_terms);
        }
    }
}

void InitializeMatrix(std::vector<double>& matrix, double p) {
    const int matrix_size = sqrt(matrix.size());
    const int num_intervals = matrix_size + 1;
    const double coeff = double((num_intervals-1) * (num_intervals-1));

    for (int i = 0; i < matrix_size; ++i) {
        for (int j = 0; j < matrix_size; ++j) {
            if (i == j) {
                matrix[i * matrix_size + j] = p + 2.0 * coeff;
            }
            else if (abs(i - j) == 1) {
                matrix[i * matrix_size + j] = -coeff;
            }
            else {
                matrix[i * matrix_size + j] = 0.0;
            }
        }
    }

    matrix[0] = 3.0 * coeff + p;
    matrix[num_intervals * (num_intervals - 2)] = 3.0 * coeff + p;
}

std::vector<double> MatrixVectorProduct(const std::vector<double>& mat,
                                      const std::vector<double>& vec) {
    const int sizee = static_cast<int>(sqrt(mat.size()));
    std::vector<double> product(sizee, 0.0);

    for (int i = 0; i < sizee; ++i) {
        for (int k = 0; k < sizee; ++k) {
            product[i] += mat[i*sizee + k] * vec[k];
        }
    }
    return product;
}

double ComputeResidualNorm(const std::vector<double>& matrix,
                         const std::vector<double>& right_hand_side,
                         const std::vector<double>& solution_vector_x,
                         int size,
                         std::vector<double>& tmp_vector) {
    double residual_squared = 0.0;

    double h = 1.0 / (double)(size - 1);
    tmp_vector = MatrixVectorProduct(matrix, solution_vector_x);

    for (int i = 0; i < size - 1; ++i) {
        double diff = right_hand_side[i] - tmp_vector[i];
        residual_squared += diff * diff;
    }

    return sqrt(residual_squared * h);
}

double FillTridiagonalMatrixWithVariableParameter(std::vector<double>& matrix) {
    int size = sqrt(matrix.size()) + 1;
    double average_value = 0.0;
    double param_k; //pk

    for (int i = 0; i < size - 1; ++i) {
        param_k = 1.0;
        for (int j = 0; j < size - 1; ++j) {
            if (i == j) {
                matrix[i * (size - 1) + j] =
                    2.0 * (double)((size - 1) * (size - 1)) + 1.0 + param_k * param_k;
            }
            else if (abs(i - j) == 1) {
                matrix[i * (size - 1) + j] =
                    -(double)((size - 1) * (size - 1));
            }
            else {
                matrix[i * (size - 1) + j] = 0.0;
            }
        }
        average_value += 1.0 + param_k * param_k;
    }

    average_value /= (double)(size - 1);
    param_k = 1.0;
    matrix[0] = 3.0 * (double)((size - 1) * (size - 1)) + param_k;
    //matrix[size * (size - 2)] = 3.0 * (double)((size - 1) * (size - 1)) + param_k;

    return average_value;
}
