#include "system.hpp"

// Metod Gaussa s vyborom glavnogo elementa po stolbtsam
int solvingSLE(std::vector<double>& A, std::vector<double>& b, std::vector<double>& Coef, std::vector<int>& memory, int N){
    for (int j = 0; j < N; j++) {
        memory[j] = j;
    }

    double mean = 1e-4;
    for (int i = 0; i < N * N; i++) {
        mean += std::fabs(A[i]);
    }

    double eps = mean * 1e-10 / (N * N);

    for (int step = 1; step <= N; step++) {
        int indMax = step;
        double max = A[e(step, step, N)]; // Diagonalnyy element

        // Ishchem maksimalnyy element po stolbtsu
        for (int i = indMax + 1; i <= N; i++) {
            if (std::fabs(A[e(step, i, N)]) > std::fabs(max)){
                indMax = i;
                max = A[e(step, i, N)];
            }
        }

        if (std::fabs(max) < std::fabs(eps)) {
            return -1;
        }

        // Menyayem stolbtsy mestami
        if (indMax != step) {
            for (int i = 1; i <= N; i++) {
                std::swap(A[e(i, step, N)], A[e(i, indMax, N)]);
            }
            std::swap(memory[indMax - 1], memory[step - 1]);
        }

        // Pryamoy khod metoda Gaussa. ustranyaya elementy pod diagonalyu
        for (int i = step + 1; i <= N; i++) {
            double v = (A[e(i, step, N)] / max);
            for (int j = step; j <= N; j++) {
                A[e(i, j, N)] -= v * A[e(step, j, N)];
            }
            b[i - 1] -= v * b[step - 1];
        }
    }

    // Obratnyy khod metoda Gaussa
    for (int step = 1; step <= N; step++) {
        for (int i = 1; i <= N - step; i++) {
            b[i - 1] -= b[N - step] * A[e(i, N - step + 1, N)] / A[e(N - step + 1, N - step + 1, N)];
            A[e(i, N - step + 1, N)] = 0;
        }
        b[N - step] /= A[e(N - step + 1, N - step + 1, N)];
        A[e(N - step + 1, N - step + 1, N)] = 1;
    }

    // Vosstanovleniye iskhodnogo poryadka peremennykh
    for (int k = 0; k < N; k++) {
        Coef[memory[k]] = b[k];
    }
    return 0;
}
