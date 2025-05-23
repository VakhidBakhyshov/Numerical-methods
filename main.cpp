#include "task_3.hpp"
#include <fstream>

int main(int argc, char *argv[]) {

    //N - number of time grid points
    //M - number of spatial grid points
    //type_of_scheme = 1 - scheme is explicit or type_of_scheme = 2 - scheme is implicit
    //key = 1 or 2 - if key = 1 - calculating convergence rate for explicit scheme M, N -> 2M, 2N)
    //if key = 2 - calculating convergence rate for implicit scheme M, N -> 2M, 4N)
    //if key = 3 - calculating convergence rate for implicit scheme M = N
    int N, M, type_of_scheme, key;

    if (argc != 5) {std::cerr<<"Error input. Please enter argc = 5!\n"; return -1;}

    if ((sscanf(argv[1], "%d", &N) != 1) || (sscanf(argv[2], "%d", &M) != 1) || (sscanf(argv[3], "%d", &type_of_scheme) != 1)
                || (sscanf(argv[4], "%d", &key)!=1)){
        std::cout<<"Invalid input!\n";
        return -1;
    }

    if (((type_of_scheme != 1) && (type_of_scheme != 2)) || (N < 3) || (M < 3) || ((key != 1) && (key != 2) && (key != 3))) {
        throw std::invalid_argument("x_grid and t_grid sizes must be >= 3; or wrong selected type of scheme");
    }

    double T = 1.0;
    double h = 1.0 / (double)(M - 1);
    double tau = T/N;

    std::vector<double> x_array(M + 1);

    Filling_x_Grid(x_array, h);

    std::ofstream file_1("x.txt"); //x_grid
    if (!file_1.is_open()) throw std::runtime_error("Failed to open file x.txt\n");

    for (int i = 0; i <= M; ++i) {
        file_1 << std::setprecision(6) << x_array[i] << "\n";
    }

    std::vector<double> t_array(N + 1);

    Filling_t_Grid(t_array, tau);

    std::ofstream file_2("t.txt"); //t_grid
    if (!file_2.is_open()) throw std::runtime_error("Failed to open file t.txt\n");

    for (int i = 0; i <= N; ++i) {
        file_2 << std::setprecision(6) << t_array[i] << "\n";
    }

    std::vector<double> u((N+1)*(M+1));

    if (type_of_scheme == 1) ExplicitScheme(N, M, tau, h, u);
    else if (type_of_scheme == 2) ImplicitScheme(N, M, tau, h, u);
    else throw std::invalid_argument("Invalid scheme, use 'explicit' or 'implicit'");

    std::ofstream file_3("matrix.txt"); //it's our matrix solution
    if (!file_3.is_open()) throw std::runtime_error("Failed to open file matrix.txt");

    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= M; ++j){
            file_3 << std::setprecision(6) << u[i*(M+1)+j];
            if (j < M) file_3 << " ";
        }
        file_3 << "\n";
    }

    std::cout << "For N = " << N << " and M = " << M << ": " << "L_2-norma = " << norma(u, N, M, h, T) << std::endl;

    if (type_of_scheme == 1 || type_of_scheme == 2) convergence_rate(T, type_of_scheme, key);

    file_1.close();
    file_2.close();
    file_3.close();

    return 0;
}