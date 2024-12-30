#include "solution.hpp"

int main(int argc, char* argv[]){

    int numberOfTests = 19;
    int typeOfQuadrature, N;
    double xa, xb, ya, yb;
    //int typeOfQuadrature = 0;
    // double sum = 0;
    if (argc<7 || argc>7){ std::cout<<"Please enter argc = 7!\n"; return -1;}
    if ((sscanf(argv[1], "%d", &N) != 1) || (N<1) || (sscanf(argv[2], "%lf", &xa)!=1) || (sscanf(argv[3], "%lf", &xb)!=1)
    || (sscanf(argv[4], "%lf", &ya)!=1) || (sscanf(argv[5], "%lf", &yb)!=1) || (sscanf(argv[6], "%d", &typeOfQuadrature)!=1)){
        std::cout<<"Invalid input!\n"<<std::endl;
        return -1;
    }
    double *VertexCoords = new double[6];
    double *CoefOfTransform = new double[6];
    double *CoordsAndWeights = new double[7 * 3];


    readingQuadratureFile(CoordsAndWeights);

    // for (int i = 0; i < 7; i++){
    //     std::cout << i+1 << " " << std::setprecision(17) << CoordsAndWeights[3*i] << " " << CoordsAndWeights[3*i + 1] << " " << CoordsAndWeights[3*i + 2] << "\n";
    // }
    std::cout << "x^k + y^k" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    //for (int k = 0; k < 50; k++){
    //int k = 13;

    //triangulationFixed(xa, xb, ya, yb, N);
    triangulation(xa, xb, ya, yb, N);

    //trueValue(xa, xb, ya, yb, k);

    for (int q = 0; q <= 33; q++) {
        for (int p = 0; p <= 33-q; p++) {
            triangulationFixed(xa, xb, ya, yb, N);
            triangulation(xa, xb, ya, yb, N);

            double true_value = trueValue(xa, xb, ya, yb, p, q);
            double approximate_value = integrateAccordingToTheQuadratureFromTheFile(f, xa, xb, ya, yb, N, VertexCoords, CoefOfTransform, CoordsAndWeights, p, q);
            double absoluteError = fabs(true_value - approximate_value);
            double relativeError = absoluteError / true_value;

            std::cout << std::left << std::setw(5) << p
                      << std::setw(5) << q
                      << std::setw(25) << std::setprecision(15) << true_value
                      << std::setw(25) << std::setprecision(15) << approximate_value
                      << std::setw(25) << std::setprecision(15) << absoluteError
                      << std::setw(25) << std::setprecision(15) << relativeError
                      << std::endl;

            writeToFileForP(f, xa, xb, ya, yb, N, numberOfTests, VertexCoords, CoefOfTransform, CoordsAndWeights, typeOfQuadrature, p, q);
        }
    }

    // case when we integrate by using file, when typeOfQuadrature = 0
    /*std::cout << "k = " << k << std::endl;

    std::cout << "True value = " << std::setprecision(15) << trueValue(xa, xb, ya, yb, k) << std::endl;

    std::cout << "Approximate value = " << std::setprecision(15)
    << integrateAccordingToTheQuadratureFromTheFile(f, xa, xb, ya, yb, k, N, VertexCoords, CoefOfTransform, CoordsAndWeights) << std::endl;

    std::cout << "Absolute error = " << std::setprecision(15)
    << fabs(trueValue(xa, xb, ya, yb, k) - integrateAccordingToTheQuadratureFromTheFile(f, xa, xb, ya, yb, k, N, VertexCoords, CoefOfTransform, CoordsAndWeights))
    << std::endl;

    std::cout << "Relative error = "
    << fabs(trueValue(xa, xb, ya, yb, k) - integrateAccordingToTheQuadratureFromTheFile(f, xa, xb, ya, yb, k, N, VertexCoords, CoefOfTransform, CoordsAndWeights)) / trueValue(xa, xb, ya, yb, k)
    << std::endl;*/



    //case when we integrate without file, when typeOfQuadrature = 1
    /*std::cout << "k = " << k << std::endl;

    std::cout << "True value = " << std::setprecision(15) << trueValue(xa, xb, ya, yb, k) << std::endl;

    std::cout << "Approximate value = " << std::setprecision(15)
    << integrateAccordingToTheQuadratureOfTheCenterOfMass(f, xa, xb, ya, yb, k, N, VertexCoords) << std::endl;

    std::cout << "Absolute error = " << std::setprecision(15)
    << fabs(trueValue(xa, xb, ya, yb, k) - integrateAccordingToTheQuadratureOfTheCenterOfMass(f, xa, xb, ya, yb, k, N, VertexCoords))
    << std::endl;

    std::cout << "Relative error = "
    << fabs(trueValue(xa, xb, ya, yb, k) - integrateAccordingToTheQuadratureOfTheCenterOfMass(f, xa, xb, ya, yb, k, N, VertexCoords)) / trueValue(xa, xb, ya, yb, k)
    << std::endl;

    writeToFileForP(f, xa, xb, ya, yb, k, N, numberOfTests, VertexCoords, CoefOfTransform, CoordsAndWeights, typeOfQuadrature);

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;*/

    //}

    delete[] VertexCoords;
    delete[] CoefOfTransform;
    delete[] CoordsAndWeights;
    return 0;
}
