#include "task_3.hpp"

int main(int argc, char* argv[]){
    std::cout << "argc = " << argc << std::endl;
    int N, k;

    //if (argc<3 || argc>5){ printf("Please enter argc 4 or 5!\n"); return -1;}
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<2) || (sscanf(argv[2], "%d", &k)!=1)){
        std::cout<<"Invalid input!\n"<<std::endl;
        return -1;
    }

    double point = 2.2;

    std::vector<double> M(N*N);
    std::vector<double> B(N) ;
    std::vector<int> memory(N);
    std::vector<double> coeffs(N);
    std::vector<double> nodes(N);
    std::vector<double> values(N);

    double a = -1., b= 1.;
    //double a = 0., b = 5.;

    GenerateEquidistantNodes(a, b, f, nodes, values);
    //GenerateChebyshevNodes(a, b, f, nodes, values);
    MatrixFill(M, nodes);
    VecFill(B, values);
    // std::cout<<"----------------------------Nodes----------------------------"<<std::endl; printVector(nodes);
    // std::cout<<"----------------------------Matrix----------------------------"<<std::endl; printMatrix(M);
    // std::cout<<"----------------------------Values----------------------------"<<std::endl; printVector(B);

    if (solve(M, B, coeffs, memory) == 0) std::cout<<"Result of Pn calculation:  "<<PnCalculation(coeffs, point)<<std::endl;
    else std::cout << "Singular matrix or error in solving." << std::endl;

    /*for (int i = 0; i < N; i++){
        std::cout << "x^"<< i << " " << coeffs[i] << std::endl;
     }*/

    LnCalculation(nodes, values, point);
    //std::cout<<"Result of Ln calculation:  "<<LnCalculation(nodes, values, point)<<std::endl;

    std::string filename = argv[3];
    std::string filename2;

    if (k == 0) WriteToFilePn(a, b, filename, coeffs);
    else if (k == 1) WriteToFileLn(a, b, filename, nodes, values);
    else if (k == 2){
	WriteToFilePn(a, b, filename, coeffs);
        filename2 = argv[4];
        WriteToFileLn(a, b, filename2, nodes, values);
    }
    else std::cout << "error\n" << std::endl;

    /*switch (k) {
        case 0:
            WriteToFilePn(a, b, filename, coeffs);
            break;
        case 1:
            WriteToFileLn(a, b, filename, nodes, values);
            break;
        case 2:
            WriteToFilePn(a, b, filename, coeffs);
            filename2 = argv[4];
            WriteToFileLn(a, b, filename2, nodes, values);
            break;
        default:
            std::cerr << "no file to write :( " << std::endl;
            break;
    }*/

    return 0;
}
