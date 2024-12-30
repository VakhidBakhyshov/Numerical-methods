#include "solution.hpp"

/*double f(double x, double y, int k){
    if (k == 0) return 1 + x * 0 + y * 0;
    else if (k == 1) return x;
    else if (k == 2) return y;
    else if (k == 3) return x * x;
    else if (k == 4) return y * y;
    else if (k == 5) return x * x + y * y;
    else if (k == 6) return x * y;
    else if (k == 7) return x * x * x;
    else if (k == 8) return y * y * y;
    else if (k == 9) return pow(x, 4);
    else if (k == 10) return pow(y, 4);
    else if (k == 11) return pow(x, 5);
    else if (k == 12) return pow(y, 5);
    else if (k == 13) return sqrt(x * y);
    else return pow(x, 4) + pow(x, 2) * pow(y, 2) + pow(y, 4);
}

double trueValue(double xa, double xb, double ya, double yb, int k){

    double truevalue = 0;
    if (k == 0) truevalue = (xb - xa) * (yb - ya);
    //if (k == 0) truevalue = (xb - xa) * (yb - ya) * (xb * xb + xb * xa + xa * xa + yb * yb + yb * ya + ya * ya) / 3.;
    else if (k == 1) truevalue = (xb * xb - xa * xa) * (yb - ya) / 2.;
    else if (k == 2) truevalue = (xb - xa) * (yb * yb - ya * ya) / 2.;
    else if (k == 3) truevalue = (xb * xb * xb - xa * xa * xa) * (yb - ya) / 3.;
    else if (k == 4) truevalue = (xb - xa) * (yb * yb * yb - ya * ya * ya) / 3.;
    else if (k == 5) truevalue = (xb - xa) * (yb - ya) * (xb * xb + xb * xa + xa * xa + yb * yb + yb * ya + ya * ya) / 3.;
    else if (k == 6) truevalue = (xb * xb - xa * xa) * (yb * yb - ya * ya) / 4.;
    else if (k == 7) truevalue = (pow(xb, 4) - pow(xa, 4)) * (yb - ya) / 4.;
    else if (k == 8) truevalue = (xb - xa) * (pow(yb, 4) - pow(ya, 4)) / 4.;
    else if (k == 9) truevalue = (pow(xb, 5) - pow(xa, 5)) * (yb - ya) / 5.;
    else if (k == 10) truevalue = (xb - xa) * (pow(yb, 5) - pow(ya, 5)) / 5.;
    else if (k == 11) truevalue = (pow(xb, 6) - pow(xa, 6)) * (yb - ya) / 6.;
    else if (k == 12) truevalue = (xb - xa) * (pow(yb, 6) - pow(ya, 6)) / 6.;
    else if (k == 13) truevalue = 4 * (pow(xb, 3 / 2.) - pow(xa, 3 / 2.)) * (pow(yb, 3 / 2.) - pow(ya, 3 / 2.)) / 9.;
    else truevalue = 0.2 * (pow(xa, 5) - pow(xb, 5)) * (ya - yb) + 1./9. * (pow(xa, 3) - pow(xb, 3)) * (pow(ya, 3) - pow(yb, 3)) +  0.2 * (pow(ya, 5) - pow(yb, 5)) * (xa - xb);
    //truevalue = ((pow(xb, k + 1) - pow(xa, k + 1)) * (yb - ya) + (pow(yb, k + 1) - pow(ya, k + 1)) * (xb - xa)) / (double)(k + 1);
    return truevalue;
}*/

double f(double x, double y, int p, int q){
    return pow(x, p) * pow(y, q);
}

double trueValue(double xa, double xb, double ya, double yb, int p, int q){
    double integralX, integralY;

    if (p == -1) {
        integralX = log(xb) - log(xa); // Для a = -1
    } else {
        integralX = (pow(xb, p + 1) - pow(xa, p + 1)) / (p + 1);
    }

    if (q == -1) {
        integralY = log(yb) - log(ya); // Для b = -1
    } else {
        integralY = (pow(yb, q + 1) - pow(ya, q + 1)) / (q + 1);
    }

    return integralX * integralY;
}

void triangulationFixed(double xa, double xb, double ya, double yb, int N){
    std::ofstream outFile("outstart.txt");

    double hx = (xb -xa) / (double)N;
    double hy = (yb -ya) / (double)N;

    outFile << "Number of nodes: " << (N + 1) * (N + 1) << std::setw(10) <<    " cord x " << std::setw(21) << "cord y " << "\n";

    for (int i = 0; i < (N + 1) * (N + 1); i++){
        outFile << std::setw(4) << i + 1 << ": "
            << std::setw(20) << xa + (i % (N + 1)) * hx << " "
            << std::setw(20) << ya + (i / (N + 1)) * hy << "\n";
    }

    double vert_left_down, vert_left_up, vert_right_down, vert_right_up;
    int triangle_index = 1;

    outFile << "Number of triangles: " << 2 * N * N << std::setw(7) << "vert1 "<< std::setw(20) << "vert2 "<< std::setw(21) << "vert3 " << "\n";

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            vert_left_down = i * (N + 1 ) + j + 1;
            vert_left_up = (i + 1) * (N + 1) + j + 1;
            vert_right_down = i * (N + 1) + (j + 1) + 1;
            vert_right_up = (i + 1) * (N + 1) + (j + 1) + 1;

            outFile << std::setw(4) << triangle_index << " "
                << std::setw(20) << vert_left_down << " "
                << std::setw(20) << vert_right_down << " "
                << std::setw(20) << vert_left_up << "\n";

            triangle_index++;  //triangle_index — schetchik treugolnikov.

            outFile << std::setw(4) << triangle_index << " "
                << std::setw(20) << vert_right_up << " "
                << std::setw(20) << vert_right_down << " "
                << std::setw(20) << vert_left_up << "\n";

            triangle_index++;
        }
    }
    outFile.close();
}

void triangulation(double xa, double xb, double ya, double yb, int N){
    std::ofstream outFile("out.txt");

    double hx = (xb -xa) / (double)N;
    double hy = (yb -ya) / (double)N;

    outFile << "Number of nodes: " << (N + 1) * (N + 1) << std::setw(10) <<    " cord x " << std::setw(21) << "cord y " << "\n";

    for (int i = 0; i < (N + 1) * (N + 1); i++){
        outFile << std::setw(4) << i + 1 << ": "
            << std::setw(20) << xa + (i % (N + 1)) * hx << " "
            << std::setw(20) << ya + (i / (N + 1)) * hy << "\n";
    }

    double vert_left_down, vert_left_up, vert_right_down, vert_right_up;
    int triangle_index = 1;

    outFile << "Number of triangles: " << 2 * N * N << std::setw(7) << "vert1 "<< std::setw(20) << "vert2 "<< std::setw(21) << "vert3 " << "\n";

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            vert_left_down = i * (N + 1 ) + j + 1;
            vert_left_up = (i + 1) * (N + 1) + j + 1;
            vert_right_down = i * (N + 1) + (j + 1) + 1;
            vert_right_up = (i + 1) * (N + 1) + (j + 1) + 1;

            outFile << std::setw(4) << triangle_index << " "
                << std::setw(20) << vert_left_down << " "
                << std::setw(20) << vert_right_down << " "
                << std::setw(20) << vert_left_up << "\n";

            triangle_index++;

            outFile << std::setw(4) << triangle_index << " "
                << std::setw(20) << vert_right_up << " "
                << std::setw(20) << vert_right_down << " "
                << std::setw(20) << vert_left_up << "\n";

            triangle_index++;
        }
    }
    outFile.close();
}

//Eta funktsiya vychislyayet koordinaty vershiny treugolnika po eye indeksu m v setke.
void getVertexCoordinates(double xa, double xb, double ya, double yb, int m, int N, double &xm, double &ym){
    //m — indeks vershiny; xm. ym — koordinaty vershiny.
    //Indeksy i, j nakhodyatsya iz m. posle chego koordinaty vychislyayutsya s uchetom shaga.
    int i = (m - 1) / (N + 1);
    int j = (m - 1) % (N + 1);

    double hx = (xb - xa) / N;
    double hy = (yb - ya) / N;

    xm = xa + j * hx;
    ym = ya + i * hy;
}

//funktsiya vypolnyayet chislennoye integrirovaniye funktsii s ispolzovaniyem kvadratury tsentra mass dlya kazhdogo treugolnika.

//double integrateAccordingToTheQuadratureOfTheCenterOfMass(double (*f)(double, double, int), double xa, double xb, double ya, double yb, int k, int N, double *VertexCoords){
double integrateAccordingToTheQuadratureOfTheCenterOfMass(double (*f)(double, double, int, int),
                                                            double xa, double xb, double ya, double yb,
                                                             int N, double *VertexCoords, int p, int q){
    //VertexCoords — massiv dlya khraneniya koordinat vershin treugolnikov.
    //res — nakoplennoye znacheniye integrala.
    double res = 0.;
    double hx = (xb-xa) / (double)N;
    double hy = (yb-ya) / (double)N;

    //Chteniye treugolnikov iz fayla out.txt
    std::ifstream file("out.txt");

    std::string line; // A variable for storing a string that we read from a file
    bool readingTriangles = false;

    //vychisleniye ikh tsentra mass i dobavleniye vklada v integral.
    while (getline(file, line)){ // Reads one line from the file and saves it to the line variable, works as long as there are lines in the file
        // нам нужна инфа после слов в файле "Number of triangles:"
        if (line.find("Number of triangles:") != std::string::npos){
            readingTriangles = true;
            continue;
        }

        if (readingTriangles){
            std::stringstream ss(line);
            int n, v1, v2, v3;
            if (ss >> n >> v1 >> v2 >> v3){

                // Calculating the coordinates of the vertices of the triangle
                getVertexCoordinates(xa, xb, ya, yb, v1, N, VertexCoords[0], VertexCoords[1]);
                // std::cout << "v1 = " << v1 << " VertexCoords[0] = " << VertexCoords[0] << " VertexCoords[1] = " << VertexCoords[1] << "\n";
                getVertexCoordinates(xa, xb, ya, yb, v2, N, VertexCoords[2], VertexCoords[3]);
                // std::cout << "v2 = " << v2 << " VertexCoords[2] = " << VertexCoords[2] << " VertexCoords[3] = " << VertexCoords[3] << "\n";
                getVertexCoordinates(xa, xb, ya, yb, v3, N, VertexCoords[4], VertexCoords[5]);
                // std::cout << "v3 = " << v3 << " VertexCoords[4] = " << VertexCoords[4] << " VertexCoords[5] = " << VertexCoords[5] << "\n";

                // The center of mass of the triangle; xm. ym — koordinaty tsentra mass.
                double xm = (VertexCoords[0] + VertexCoords[2] + VertexCoords[4]) / 3.0;
                double ym = (VertexCoords[1] + VertexCoords[3] + VertexCoords[5]) / 3.0;

                // Calculating: the area of the triangle * the value of the function
                //res += 0.5 * hx * hy * f(xm, ym, k);
                res += 0.5 * hx * hy * f(xm, ym, p, q);
            }
        }
    }
    file.close();
    return res;
}

void readingQuadratureFile(double *CoordsAndWeights){
    //std::ifstream file("CoordsAndWeights.txt");
    std::ifstream file("Coords_And_Weights.txt");

    int i = 0;
    std::string line; // A variable for storing a string that we read from a file
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double x, y, weight;
        while (iss >> x >> y >> weight) {
            CoordsAndWeights[3*i] = x;
            CoordsAndWeights[3*i+1] = y;
            CoordsAndWeights[3*i+2] = weight;
            i++;
        }
    }
    file.close();
}

void searchAffineTransformation(double* VertexCoords, double *EquationCoefficients){
    EquationCoefficients[0] = VertexCoords[2] - VertexCoords[0]; // a
    EquationCoefficients[1] = VertexCoords[4] - VertexCoords[0]; // b
    EquationCoefficients[2] = VertexCoords[3] - VertexCoords[1]; // c
    EquationCoefficients[3] = VertexCoords[5] - VertexCoords[1]; // d
    EquationCoefficients[4] = VertexCoords[0]; // e
    EquationCoefficients[5] = VertexCoords[1]; // f
}

//double integrateAccordingToTheQuadratureFromTheFile(double (*f)(double, double, int), double xa, double xb, double ya, double yb, int k, int N, double *VertexCoords, double *CoefOfTransform, double *CoordsAndWeights){
double integrateAccordingToTheQuadratureFromTheFile(double (*f)(double, double, int, int), double xa, double xb, double ya, double yb, int N, double *VertexCoords, double *CoefOfTransform, double *CoordsAndWeights, int p, int q){

    double res = 0.;
    double hx = (xb-xa) / (double)N;
    double hy = (yb-ya) / (double)N;
    double x, y;

    std::ifstream file("out.txt");

    std::string line; // A variable for storing a string that we read from a file
    bool readingTriangles = false;
    int t = 0;

    while (getline(file, line)){ // Reads one line from the file and saves it to the line variable, works as long as there are lines in the file
        if (line.find("Number of triangles:") != std::string::npos){
            readingTriangles = true;
            continue;
        }

        if (readingTriangles){
            std::stringstream ss(line);
            int n, v1, v2, v3;
            if (ss >> n >> v1 >> v2 >> v3){

                // Calculating the coordinates of the vertices of the triangle
                getVertexCoordinates(xa, xb, ya, yb, v1, N, VertexCoords[0], VertexCoords[1]);
                getVertexCoordinates(xa, xb, ya, yb, v2, N, VertexCoords[2], VertexCoords[3]);
                getVertexCoordinates(xa, xb, ya, yb, v3, N, VertexCoords[4], VertexCoords[5]);
                // if (t < 2){
                //     std::cout << "t = " << t << "\n";
                //     std::cout << "v1 = " << v1 << " VertexCoords[0] = " << VertexCoords[0] << " VertexCoords[1] = " << VertexCoords[1] << "\n";
                //     std::cout << "v2 = " << v2 << " VertexCoords[2] = " << VertexCoords[2] << " VertexCoords[3] = " << VertexCoords[3] << "\n";
                //     std::cout << "v3 = " << v3 << " VertexCoords[4] = " << VertexCoords[4] << " VertexCoords[5] = " << VertexCoords[5] << "\n";
                // }

                searchAffineTransformation(VertexCoords, CoefOfTransform);

                for (int i = 0; i < 7; i++){
                    // if (t == 0 && i == 0) std::cout << "i = " << i << " a11 = " << CoefOfTransform[0] << " a12 = " << CoefOfTransform[1] << " a21 = " << CoefOfTransform[2] << " a22 = " << CoefOfTransform[3] << " b1 = " << CoefOfTransform[4] << " b2 = " << CoefOfTransform[5] << "\n";
                    x = CoefOfTransform[0] * CoordsAndWeights[3*i] + CoefOfTransform[1] * CoordsAndWeights[3*i + 1] + CoefOfTransform[4]; // x = a * x + b * y + e; x = CoordsAndWeights[3*i];
                    y = CoefOfTransform[2] * CoordsAndWeights[3*i] + CoefOfTransform[3] * CoordsAndWeights[3*i + 1] + CoefOfTransform[5]; // y = c * x + d * y + f; y = CoordsAndWeights[3*i + 1];
                    // if (t == 0) std::cout << std::left  << "x = " << std::fixed << std::setw(19) << std::setprecision(17) << x << " y = " << y << "\n";
                    //res += CoordsAndWeights[3*i + 2] * f(x, y, k); // (CoefOfTransform[0] * CoefOfTransform[3] - CoefOfTransform[1] * CoefOfTransform[2]);
                    res += CoordsAndWeights[3*i + 2] * f(x, y, p, q);
                }
            }
            t++;
        }
    }
    res = res * hx * hy * 0.5 * 2;
    file.close();
    return res;
}

//void writeToFileForP(double (*f)(double, double, int), double xa, double xb, double ya, double yb, int k, int N, int numTests, double *VertexCoords, double *CoefOfTransform, double *CoordsAndWeights, int typeOfQuadrature){
void writeToFileForP(double (*f)(double, double, int, int),
                        double xa, double xb, double ya, double yb, int N,
                        int numTests, double *VertexCoords, double *CoefOfTransform,
                        double *CoordsAndWeights, int typeOfQuadrature, int p, int q){
    double approximate_value, true_value, error;

    std::ofstream outFile("p.txt");
    outFile<<std::setw(10) << std::left
         << std::setw(15) << "N"
         << std::setw(20) << "Aproximate"
         << std::setw(20) << "True"
         << std::setw(20) << "Error"
         << std::endl;

    for(int i = 0; i < numTests; i++){
        triangulation(xa, xb, ya, yb, N);
        if (typeOfQuadrature == 0){
            //approximate_value = integrateAccordingToTheQuadratureOfTheCenterOfMass(f, xa, xb, ya, yb, k, N, VertexCoords);
            approximate_value = integrateAccordingToTheQuadratureOfTheCenterOfMass(f, xa, xb, ya, yb, N, VertexCoords, p, q);
        }
        else if (typeOfQuadrature == 1){
            //approximate_value = integrateAccordingToTheQuadratureFromTheFile(f, xa, xb, ya, yb, k, N, VertexCoords, CoefOfTransform, CoordsAndWeights);
            approximate_value = integrateAccordingToTheQuadratureFromTheFile(f, xa, xb, ya, yb, N, VertexCoords, CoefOfTransform, CoordsAndWeights, p, q);
        }
        else{
            std::cout << "Type of quadrature is not correct" << std::endl;
            break;
        }
        //true_value = trueValue(xa, xb, ya, yb, k);
        true_value = trueValue(xa, xb, ya, yb, p, q);
        error = fabs(approximate_value - true_value);

            outFile << std::setw(10) << std::left
                << std::setw(6) << std::fixed << std::setprecision(3) << N
                << std::setw(24) << std::fixed << std::setprecision(20) << approximate_value
                << std::setw(24) << std::fixed << std::setprecision(20) << true_value
                << std::setw(24) << std::fixed << std::setprecision(20) << error
                << std::endl;

            N++;
    }
    outFile.close();
}
