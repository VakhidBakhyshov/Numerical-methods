#include "task_2.hpp"


int main(int argc, char *argv[]){

    int N, k;

    //if (argc<3 || argc>4){ printf("Please enter argc 3 or 4!\n"); return -1;}
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<3) || (sscanf(argv[2], "%d", &k)!=1) || ((k==0) && (argc==3))){
        std::cout<<"Invalid input! \n"<<std::endl;
        return -1;
    }

    double* xk      = new double[N+1];
    double* fmemory = new double[N+1];
    double* U       = new double[(N+1)*(N+1)];
    double* D       = new double[(N+1)*(N+1)];
    double* C       = new double[(N+1)*(N+1)];
    double* phi     = new double[N+1];

    /*if ((!xk) || (!U) || (!D) || (!C) || (!phi) || (!fmemory)){
        std::cout<<"Not enough memory!\n";
                return -1;
    }*/

    if (k > 0){
        std::string filename = argv[3];
        WriteToFile(filename, N, xk, U, C, D, fmemory, phi);
        pcalculate(N);
    }
    else{
        std::cout<<WriteToConsole(N, xk, U, C, D, fmemory, phi)<< " = DURATION (fToC)\n" << std::endl;
    }

    delete[] xk;
    delete[] fmemory;
    delete[] U;
    delete[] D;
    delete[] C;
    delete[] phi;

    return 0;
}