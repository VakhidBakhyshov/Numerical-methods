#include "task_4.hpp"

double f(double x){
    //return((exp(x)-exp(1.0))*sin(x));
    //return ((x > 0.3) && (x < 0.5) ? 1. : 0.);
    //return exp(x);
    //return 16*pow(x,5);
    //return (x*x+sin(x))*cos(3*x);
    //return 1./(1. + 25.*x*x); //Runge's function
    return std::fabs(x);
}

void fillingNodes(double a, double b, std::vector<double>& Nodes, int type){
    int N = Nodes.size();
    double now = a;
    //Generate Equidistant Nodes
    if (type == 0){
        double delta = (b - a) / (N - 1);
        for (int i = 0; i < N; i++){
	    Nodes[i] = now;
	    now += delta;
            //Nodes[i] = a + i * delta;
        }
    }
    //Generate Chebyshev Nodes
    else if (type == 1){
        double delta1 = (b + a) / 2;
        double delta2 = (b - a) / 2;
        for (int i = 0; i < N; i++){
            Nodes[N - 1 - i] = delta1 + delta2 * cos((2 * i + 1) * M_PI / (2 * N));
        }
    }
    //Generate Random nodes
    else{
        for (int i = 0; i < N; i++){
            Nodes[i] = 0;
            //srand(1);

            int randint;
            for(int i = 0; i < N; i++){
                randint = rand();
                Nodes[i] = a + (b - a) * ((double) randint / RAND_MAX);
            }
            bubbleSort(Nodes, N);
            }
    }
}

void bubbleSort(std::vector<double>& Nodes, int N) {
    for (int i = 0; i < N - 1; i++) {
        for (int j = 0; j < N - 1; j++) {
            if (Nodes[j] > Nodes[j + 1]) {
                std::swap(Nodes[j], Nodes[j + 1]);
            }
        }
    }
}

void fillingValuesInNodes(std::vector<double>& Nodes, std::vector<double>& ValuesInNodes, int N){
    for (int i = 0; i < N; i++){
    	ValuesInNodes[i] = f(Nodes[i]);
    }
}

void fillingMatrixA(std::vector<double>& A, std::vector<double>& Nodes, int N){
    for (int i = 0; i < N; i++){
        for (int j = 1; j < N; j++){
       		A[i * N + j] = pow(Nodes[i], j-1);
	}
        A[i * N] = (i % 2 == 0) ? 1.0 : -1.0;
    }
}

void printMatrixAToConsole(std::vector<double>& A, int N){
    std::cout << "========================================================================= " << " Matrix A " << " =========================================================================" << std::endl;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            std::cout << std::setw(12) << std::setprecision(4) << A[i * N + j] << " ";
        }
        std::cout << " " << std::endl;
    }
    std::cout << " " << std::endl;
}

void fillingVectorB(std::vector<double>& B, std::vector<double>& ValuesInNodes, int N){
    for (int i = 0; i < N; i++){
        B[i] = ValuesInNodes[i];
    }
}

void printCoefToConsole(std::vector<double>& Coef, int N){
    for (int i = 0; i < N; i++){
        std::cout << Coef[i] << " ";
    }
    std::cout << std::endl;
}

double pnCalculating(std::vector<double>& Coef, int N, double x){
    double res = 0;
    for (int i = 1; i < N; i++){
        res += Coef[i] * pow(x, i-1);
    }
    return res;
}

/*double CalcPolynom(std::vector<double>& Coef, double x, int N){
    double ans = coeffs[1];
    for (int i = 1; i < N - 1; ++i){
        double xi = pow(x, i);
        ans += xi * coeffs[i + 1];
    }
    return ans;
}*/

void shreddedGrid(std::vector<double>& ShereddedGrid, double *Nodes, int N){
    double h;
    for (int i = 0; i < N; i++){
        ShereddedGrid[3*i] = Nodes[i];
    }
    for(int i = 0; i < N-1; i++){
        h = (Nodes[i+1]-Nodes[i])/3;
        ShereddedGrid[3*i + 1] = ShereddedGrid[3*i] + h;
        ShereddedGrid[3*i + 2] = ShereddedGrid[3*i + 1] + h;
    }
}

void fillingShreddedGrid(std::vector<double>& ShreddedGrid, std::vector<double>& Nodes, int N){
    double h;
    for (int i = 0; i < N; i++){
        ShreddedGrid[3*i] = Nodes[i];
    }
    for(int i = 0; i < N-1; i++){
        h = (Nodes[i+1]-Nodes[i])/3;
        ShreddedGrid[3*i + 1] = ShreddedGrid[3*i] + h;
        ShreddedGrid[3*i + 2] = ShreddedGrid[3*i + 1] + h;
    }
}

void fillingValuesInShreddedGrid(std::vector<double>& ValuesInShreddedGrid, std::vector<double>& ShreddedGrid, int N){
    for(int i = 0; i < 3*N - 2; i++){
        ValuesInShreddedGrid[i]=f(ShreddedGrid[i]);
    }
}

void writeToFilePn(std::vector<double>& Coef, int N, std::vector<double>& ShreddedGrid, int altcount){
    std::ofstream outFile("Pn.txt");
    double trueValue, approximateValue, err;

    outFile<<std::setw(10)<<" "<<"x"<<std::setw(10)<<" "\
    <<std::setw(9)<<" "<<"f(x)"<<std::setw(9)<<" "\
    <<std::setw(9)<<" "<<"Pn"<<std::setw(9)<<" "\
    <<std::setw(9)<<" "<<"err"<<std::setw(9)<<" "<<std::endl;

    for (int i = 0; i < 3*N - 2; i++){
        trueValue = f(ShreddedGrid[i]);
        approximateValue =  pnCalculating(Coef, altcount, ShreddedGrid[i]);
        err = std::fabs(trueValue-approximateValue);
        outFile << std::setprecision(15) << std::fixed \
        << std::setw(20) << ShreddedGrid[i] << " " \
        << std::setw(20) << trueValue << " " \
        << std::setw(20) << approximateValue << " " \
        << std::setw(20) << err  << std::endl;
    }
    outFile.close();
}

void fillingSigma(std::vector<double>& Sigma, std::vector<double>& Nodes, int altcount, int N){
    int k;
    k = N / altcount;
    for(int i = 0; i < altcount; i++){
        Sigma[i] = Nodes[i*k];
    }
}

void fillingValuesInSigma(std::vector<double>& sigma, std::vector<double>& ValuesInSigma, int altcount){
    for (int i = 0; i < altcount; i++){
    	ValuesInSigma[i] = f(sigma[i]);
    }
}

bool maxDeviation(std::vector<double>& Nodes, std::vector<double>& sigma, std::vector<double>& Coef, std::vector<double>& ValuesInNodes, std::vector<double>& ValuesInSigma, int N, int altcount){
    double h = Coef[0]; // h - tekushchaya velichina maksimalnogo otkloneniya iz koeffitsiyentov
    // std::cout << "h = " << h << std::endl;

    double maxx = -1.; //maxx - ispolzuyetsya dlya khraneniya tekushchego maksimalnogo znacheniya
    double res = 0.; //res - vremennaya peremennaya dlya vychisleniya otkloneniya
    double res_old, res_new;
    int k = 0; //k - indeks uzla s maksimalnym otkloneniyem

    //We are looking for the maximum deviation in the nodes
    //poisk uzla s maksimalnym otkloneniyem
    for(int i = 0; i < N; i++){
        // otkloneniye |P_n(x_i) - y(x_i)|
        res = std::fabs(pnCalculating(Coef, altcount, Nodes[i]) - ValuesInNodes[i]);
        if (res > maxx){
            maxx = std::fabs(res);
            k = i;
        }
    }

    //proverka usloviya optimalnosti
    if (std::fabs(h) + 1e-10 > maxx){
                return 0;
    }
    else{
        if (Nodes[k] < sigma[0]){
        // std::cout << "case LEFT" << std::endl;
            res_old = pnCalculating(Coef, altcount, sigma[0]);
            res_new = pnCalculating(Coef, altcount, Nodes[k]);
            if (((res_old - ValuesInSigma[0] < 0) && (res_new - ValuesInNodes[k] < 0)) || ((res_old - ValuesInSigma[0] > 0) && (res_new - ValuesInNodes[k] > 0))){
                sigma[0] = Nodes[k];
                ValuesInSigma[0] = ValuesInNodes[k];
            }
            else{
                for (int i = 0; i < altcount - 1; i++){
                    sigma[altcount - i - 1] = sigma[altcount - 2 - i];
                    ValuesInSigma[altcount - i - 1] = ValuesInSigma[altcount - 2 - i];
                }
                sigma[0] = Nodes[k];
                ValuesInSigma[0] = ValuesInNodes[k];
            }
            return 1;
        }
        else if (Nodes[k] > sigma[altcount - 1]){
        // std::cout << "case RIGHT" << std::endl;

            res_old = pnCalculating(Coef, altcount, sigma[altcount - 1]);
            res_new = pnCalculating(Coef, altcount, Nodes[k]);

            if (((res_old - ValuesInSigma[altcount-1] < 0) && (res_new - ValuesInNodes[k] < 0) ) || ( (res_old - ValuesInSigma[altcount-1] > 0) && (res_new - ValuesInNodes[k] > 0))){
                sigma[altcount - 1] = Nodes[k];
                ValuesInSigma[altcount - 1] = ValuesInNodes[k];
            }
            else{
                for(int i = 0; i < altcount - 1; i++){
                    sigma[i] = sigma[i + 1];
                    ValuesInSigma[i] = ValuesInSigma[i + 1];
                }
                sigma[altcount - 1] = Nodes[k];
                ValuesInSigma[altcount - 1] = ValuesInNodes[k];
            }
            return 1;
        }
        else{
        // std::cout << "case MIDDLE" << std::endl;

            int m = 0;
            while(sigma[m] < Nodes[k]){
		        m++;
            }

            res_old = pnCalculating(Coef, altcount, sigma[m - 1]);
            res_new = pnCalculating(Coef, altcount, Nodes[k]);

            if (((res_old - ValuesInSigma[m-1] < 0) && (res_new - ValuesInNodes[k] < 0)) || ((res_old - ValuesInSigma[m-1] > 0) && (res_new - ValuesInNodes[k] > 0))){
                sigma[m - 1] = Nodes[k];
                ValuesInSigma[m - 1] = ValuesInNodes[k];
            }
            else{
                sigma[m] = Nodes[k];
                ValuesInSigma[m] = ValuesInNodes[k];
            }
            return 1;
        }
    }
}

void printVectorToConsole(std::vector<double>& Vector, int N){
    for (int i = 0; i < N; i++){
        std::cout << Vector[i] << std::endl;
    }
    std::cout << std::endl;
}
