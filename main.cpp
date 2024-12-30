#include "task_4.hpp"
#include "system.hpp"

int main(int argc, char *argv[]){

    // N - number of grid nodes
    //deg - degree of polynomial
    // type = 0 - Equidistant, 1 - Chebyshev, 2 - Random
    int N, deg, type;
    //[a,b] - segment
    double a, b;

    /*for example: int N = 10;
    int deg = 7;
    int altcount = deg + 2;
    double a = -1.0;
    double b = 1.0;
    int type = 2;*/

    if (argc<6 || argc>7){ std::cout<<"Please enter argc 6 or 7!\n"; return -1;}
    if ((sscanf(argv[1], "%d", &N) != 1) || (sscanf(argv[2], "%d", &deg) != 1) ||
     (N<2) || (sscanf(argv[3], "%lf", &a)!=1) || (sscanf(argv[4], "%lf", &b)!=1) ||
     (sscanf(argv[5], "%d", &type)!=1)){
        std::cout<<"Invalid input!\n";
        return -1;
    }

    // Number of nodes required for the algorithm De la Vallée-Poussin polynomial of degree deg, мощность альтернанса
    int altcount = deg + 2;

    std::vector<double> Nodes(N);
    std::vector<double> ValuesInNodes(N);
    std::vector<double> A(altcount * altcount);
    std::vector<double> B(altcount);
    std::vector<double> Coef(altcount);
    std::vector<int> memory(altcount);
    std::vector<double> sigma(altcount);
    std::vector<double> ValuesInSigma(altcount);
    std::vector<double> ShreddedGrid(3*N - 2);
    std::vector<double> ValuesInShreddedGrid(3*N - 2);

    if (N == altcount){ //sluchay Chebysheva: stroitsya polinom cherez resheniye sistemu lineynykh uravneniy
    //if (deg == N - 2){

        fillingNodes(a, b, Nodes, type);
        // std::cout << "Nodes = " << std::endl;
        // printVectorToConsole(Nodes , N);

        fillingValuesInNodes(Nodes, ValuesInNodes, N);
        // std::cout << "VIN = " << std::endl;
        // printVectorToConsole(ValuesInNodes , N);

        fillingMatrixA(A, Nodes, N);
        // printMatrixAToConsole(A, N);

        fillingVectorB(B, ValuesInNodes, N);
        // std::cout << "B = " << std::endl;
        // printVectorToConsole(B , N);

        solvingSLE(A, B, Coef, memory, N);
        // std::cout << "h + Coef = " << std::endl;
        // printVectorToConsole(Coef , N);

        fillingShreddedGrid(ShreddedGrid, Nodes, N);
        // std::cout << "Shredded Grid = " << std::endl;
        // printVectorToConsole(ShreddedGrid , 3*N - 2);

        fillingValuesInShreddedGrid(ValuesInShreddedGrid, ShreddedGrid, N);
        // std::cout << "VISG = " << std::endl;
        // printVectorToConsole(ValuesInShreddedGrid , 3*N - 2);

        writeToFilePn(Coef, N, ShreddedGrid, altcount);

        std::cout << "h = " << fabs(Coef[0]) << std::endl;
    }
    else { // ispolzuyetsya algoritm Valle-Pussena dlya nakhozhdeniya ravnomernogo priblizheniya
    //else if (deg < N - 2){

        fillingNodes(a, b, Nodes, type);
        // std::cout << "Nodes = " << std::endl;
        // printVectorToConsole(Nodes , N);

        fillingValuesInSigma(Nodes, ValuesInNodes, N);

        fillingSigma(sigma, Nodes, altcount, N);
        // std::cout << "Sigma = " << std::endl;
        // printVectorToConsole(sigma , altcount);

        fillingValuesInSigma(sigma, ValuesInSigma, altcount);
        // std::cout << "VIS = " << std::endl;
        // printVectorToConsole(ValuesInSigma, altcount);

        fillingMatrixA(A, sigma, altcount);
        // printMatrixAToConsole(A, altcount);

        fillingVectorB(B, ValuesInSigma, altcount);
        // std::cout << "B = " << std::endl;
        // printVectorToConsole(B , altcount);

        solvingSLE(A, B, Coef, memory, altcount);
        // std::cout << "h + Coef = " << std::endl;
        // printVectorToConsole(Coef, altcount);

        bool flag = 1;
        int  iteration = 1;

        while((flag) && (iteration < MAX_ITERATIONS)){
	//while((flag)){
            flag = maxDeviation(Nodes, sigma, Coef, ValuesInNodes, ValuesInSigma, N, altcount);

            if (flag){
                fillingMatrixA(A, sigma, altcount);
                //printMatrixAToConsole(A, altcount);

                fillingVectorB(B, ValuesInSigma, altcount);
                solvingSLE(A, B, Coef, memory, altcount);

		/*for (int i = 0; i < N; i++){
			std::cout << std::left << std::setw(20) << std::setprecision(15) << f(Nodes[i]) - pnCalculating(Coef, altcount, Nodes[i]) << " ";
		}
		std::cout << std::endl;*/

		//printVectorToConsole(ValuesInNodes, N);
		//printMatrixAToConsole(A, altcount);

                iteration += 1;
            }
        }
	for (int i = 0; i < N; i++){
        	std::cout << std::left << std::setw(20) << std::setprecision(15) << f(Nodes[i]) - pnCalculating(Coef, altcount, Nodes[i]) << "\n";
        }
        std::cout << std::endl;

	printVectorToConsole(sigma , altcount);

        std::cout << "Number of iterations = " << iteration << std::endl;

        //printCoefToConsole(Coef, altcount);
        fillingShreddedGrid(ShreddedGrid, Nodes, N);
        // std::cout << "Shredded Grid = " << std::endl;
        // printVectorToConsole(ShreddedGrid , 3*N - 2);

        fillingValuesInShreddedGrid(ValuesInShreddedGrid, ShreddedGrid, N);
       	// std::cout << "VISG = " << std::endl;
        // printVectorToConsole(ValuesInShreddedGrid , 3*N - 2);

        writeToFilePn(Coef, N, ShreddedGrid, altcount);
        std::cout << "h = " << fabs(Coef[0]) << std::endl;

	//printVectorToConsole(ValuesInSigma, altcount);
    }

    return 0;
}
