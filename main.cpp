#include "task_1.hpp"

int main(int argc, char *argv[]){
    int N;
    double p;

    if (argc<3 || argc>3){ std::cout<<"Please enter argc = 3!\n"; return -1;}
    if ((sscanf(argv[1], "%d", &N) != 1) || (sscanf(argv[2], "%lf", &p)!=1) || (N<1) || (p<0)){
        std::cout<<"Invalid input!\n";
        return -1;
    }

    int count = 8, N_tmp = N;
    std::vector<double> norma(count);
    std::vector<double> h(count);
    std::vector<double> rate(count-1);

    std::ofstream outFile("output.txt");
    if (outFile.is_open()) {
        for (int i = 0; i < count; ++i){

                h[i] = 1.0 / (double)(N - 1);

                std::vector<double> grid_nodes(N - 1);
                std::vector<double> function_values(N - 1);

                ComputeGridPoints(grid_nodes);
                CalculateFunctionValues(grid_nodes, function_values, f);

                std::vector<double> matrix((N - 1) * (N - 1));

                InitializeMatrix(matrix, p);

                std::vector<double> solution_vector_x(N - 1);

                MethodFourier(solution_vector_x, p, function_values);

                std::vector<double> tmp(N - 1);

                norma[i] = ComputeResidualNorm(matrix, function_values, solution_vector_x, N, tmp);

                if (i > 0){
                        rate[i-1] = fabs(log(norma[i - 1] / norma[i]) / log(h[i - 1] / h[i]));
                        //std::cout << "For N = " << N << ", norma = " << norma[i] << ", convergence index = " << rate[i] <<std::endl;
                        outFile << N << " " << norma[i] << " " << rate[i-1] << std::endl;
                }
                else outFile << N << " " << norma[i] << std::endl;

                N *= 2;
        }
    }

    else {
        std::cerr << "Unable to open file: " << std::endl;
        return -1;
    }

    std::cout << "For N = " << N_tmp << ", norma = " << norma[0] << std::endl;
    N_tmp *= 2;

    for (int i = 0; i < count - 1; ++i){
        std::cout << "For N = " << N_tmp << ", norma = " << norma[i+1] << ", convergence index = " << rate[i] <<std::endl;
        N_tmp *= 2;
    }

    double result = fabs(log(norma[count-2] / norma[count-1]) / log(h[count-2] / h[count-1]));
    //double result = fabs(log(norma[0] / norma[1]) / log(h[0] / h[1]));
    std::cout << "Convergence index = " << result <<std::endl;

    outFile.close();

    return 0;
}