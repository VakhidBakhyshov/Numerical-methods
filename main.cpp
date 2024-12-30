#include "task_7.hpp"

int main(int argc, char *argv[]){
    int N;
    double p;

    if (argc<3 || argc>3){ std::cout<<"Please enter argc = 3!\n"; return -1;}
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<1) || (sscanf(argv[2], "%lf", &p)!=1) || (p<0)){
        std::cout<<"Invalid input!\n";
        return -1;
    }

    std::vector<double> Nodes(N + 1);
    std::vector<double> BasicNodes(N - 1);

    std::vector<double> ValuesInNodes(N + 1);
    std::vector<double> ValuesInBasicNodes(N - 1);

    std::vector<double> MatrixA((N - 1) * (N - 1));
    std::vector<double> b(N - 1);

    std::vector<double> BasicMatrix((N - 1) * (N - 1));
    std::vector<double> BasicMatrixA((N - 1) * (N - 1));
    std::vector<double> BasicMatrixB((N - 1) * (N - 1));
    std::vector<double> FullMatrix((N + 1) * (N + 1));
    std::vector<double> FullMatrix2((N + 1) * (N + 1));

    std::vector<double> x(N - 1);
    std::vector<double> Fullx(N + 1);
     std::vector<double> Fullx2(N + 1);

    std::vector<double> mem(N - 1);
    std::vector<double> mem1(N - 1);

    FillingNodes(Nodes);
    FillingBasicNodes(BasicNodes);

    FillingValues(Nodes, ValuesInNodes, f);
    FillingValues(BasicNodes, ValuesInBasicNodes, f);

    BasisMatrixFill(p, BasicMatrix);
    printMatrix(BasicMatrix);
    std::cout << "\n" << std::endl;

    FullMatrixFill(p, FullMatrix);
    printMatrix(FullMatrix);

    Fourier(x, p, ValuesInBasicNodes);
    Fullx[0] = Nodes[0];
    for (int i = 0; i < N-1;++i){
        Fullx[i+1] = x[i];
    }
    Fullx[N] = Nodes[N];

    /*std::cout<<"check: "<<std::endl; //<- correct
    printVector(MultiplyMatrixByVector(BasicMatrix, x));
    printVector(ValuesInBasicNodes);
    std::cout << "\n" << std::endl;*/
    printVector(MultiplyMatrixByVector(FullMatrix, Fullx));
    printVector(ValuesInNodes);
    std::cout << "\n" << std::endl;
    //printVector(x);

    double EigenValueMin = Lambda(1, N, p);
    double EigenValueMax = Lambda(N - 1, N, p);
    double tau;

    tau = 2. / (EigenValueMin + EigenValueMax);
    double q = (EigenValueMax - EigenValueMin) / (EigenValueMax + EigenValueMin);

    std::ofstream fout("out.txt");
    if (!fout) {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }

    BasisMatrixFillWithVariableP(p, BasicMatrixA); //A
    FullMatrixFillWithVariableP(FullMatrix2);

    std::vector<double> x_one(N - 1);
    for (int i = 0; i < N - 1; ++i){
        x_one[i] = 1.;
    }
    std::vector<double> x_one_full(N + 1);
    for (int i = 0; i < N + 1; ++i){
        x_one_full[i] = 1.;
    }
    /*std::cout<<"variable pk\n"<<pp<<std::endl;
    printMatrix(BasicMatrixA);
    std::cout<<"\n"<<pp<<std::endl;
    printMatrix(FullMatrix2);

    std::cout<<"check one: "<<std::endl; //<- correct

    ValuesInBasicNodes = MultiplyMatrixByVector(BasicMatrixA, x_one);
    std::cout << "vector f we know\n" << std::endl;
    printVector(ValuesInNodes);

    std::cout << "vector y we know it's vector of all 1\n" << std::endl;
    printVector(x_one);
    std::vector<double> x_new(N - 1);
    std::cout << "vector y we calculate using fourier\n" << std::endl;
    printVector(x_new);

    Fourier(x_new, p, ValuesInBasicNodes);

    Fullx2[0] = 1.;
    for (int i=0;i<N-1;++i){
        Fullx2[i+1] = x_new[i];
    }
    Fullx2[N] = 1.;

    ValuesInNodes = MultiplyMatrixByVector(FullMatrix2, Fullx2);
    std::cout << "vector f we calculate using fourier\n" << std::endl;
    printVector(ValuesInNodes);*/

    /*for (int i = 0; i < N - 1; ++i){
        ValuesInBasicNodes[i] = 0.;
    }
    ValuesInBasicNodes = MultiplyMatrixByVector(BasicMatrixA, x_one);
    std::cout << "vector f we know\n" << std::endl;
    printVector(ValuesInBasicNodes);
    std::cout << "vector y we know it's vector of all 1\n" << std::endl;
    printVector(x_one);
    std::vector<double> x_new(N - 1);
    Fourier(x_new, p, ValuesInBasicNodes);
    std::cout << "vector y we calculate using fourier\n" << std::endl;
    printVector(x_new);
    std::cout << "error between vector y we found in fourier and vector y of all 1\n" << std::endl;
    for (int k = 0; k < N - 1; ++k){
        std::cout<<"with k = "<< k << ": y[k] - 1 = " << x_new[k] - 1. <<std::endl;
    }
    for (int i = 0; i < N - 1; ++i){
        ValuesInBasicNodes[i] = 0.;
    }
    ValuesInBasicNodes = MultiplyMatrixByVector(BasicMatrixA, x_new);
    std::cout << "vector f we calculate using fourier\n" << std::endl;
    printVector(ValuesInBasicNodes);*/

    BasisMatrixFill(p, BasicMatrixB); //B

    std::cout<<"Basis Matrix Fill With Variable P\n"<<std::endl;
    printMatrix(BasicMatrixA);
    //std::cout << "\n" << std::endl;
    //printMatrix(BasicMatrixB);

    int numberTest = 10000;
    double q0 = q;
    //tau = 1.;
    //tau = 0.5;
    //tau = 0.3;

    double resid0 = BSolver(x, BasicMatrixA, BasicMatrixB, ValuesInBasicNodes, tau, 1, mem, mem1, p); //for x^0
    fout << 1 << " " << resid0 << " " << q0 * resid0 << "\n";

    /*for(int i = 0; i < N-1; i++){
        ValuesInBasicNodes[i] = 1.;
    }*/

    //fout << std::fixed << std::setprecision(15);
    for (int iter = 2; iter < numberTest+1; iter = iter * 2){
    //for (int iter = 1; iter < numberTest+1; iter += 1){
        double resid = BSolver(x, BasicMatrixA, BasicMatrixB, ValuesInBasicNodes, tau, iter, mem, mem1, p);
        fout << iter << " " << resid << " " << q0 * resid0 << "\n";
        q0 *= q0;
    }

    fout.close();
    return 0;
}
