#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "task_2.hpp"

int main(int argc, char *argv[]) {

    int N; //N = 1000
    double A; // A = 1000.0

    if (argc<3 || argc>3){ std::cout<<"Please enter argc = 3!\n"; return -1;}
    if ((sscanf(argv[1], "%d", &N) != 1) || (sscanf(argv[2], "%lf", &A)!=1)){
        std::cout<<"Invalid input!\n";
        return -1;
    }

    double error1, error2, error3, error4;
    int count_test = 5;
    double y0 = 1.0;

    std::ofstream outFile("output.txt");
    if(outFile.is_open()){
        for(int i = 0; i < count_test; i++){

            error1 = ExplicitEulerMethod(y0, A, N);
            error2 = ImplicitEulerMethod(y0, A, N);
            error3 = TrapezoidalApproach(y0, A, N);
            error4 = LeapfrogTechnique(y0, A, N);
            outFile << std::setprecision(15) << N << " " << error1 << " " << error2 << " " << error3 << " " << error4 << std::endl;
            N *= 2;

        }
    } else {
        std::cerr << "error\n" << std::endl;
        return 1;
    }
    outFile.close();

    return 0;
}
