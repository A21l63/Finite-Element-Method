#include <iostream>
#include <cmath>
#include "Equation.cpp"
#include "Solver.h"
#include "const.h"

int main() {
    Solver solver = Solver(NUM, LEFT, RIGHT);

    solver.SolveGal();
    /*solver.SolveFDM(1);


    double y_orig;
    std::cout << "x |  y[i] | sin(PI * x) | * NUM";
    for (int i = 0; i <= NUM; i++) {
        y_orig = solver.finiteDiff.grid.X[i] * (1  - solver.finiteDiff.grid.X[i]) / 2;; //X[i] * cos(M_PI * X[i]); //sin(M_PI * X[i]);
        std::cout << solver.finiteDiff.grid.X[i] << " | " << solver.finiteDiff.y[i] << " | " << y_orig << " | " << std::abs(y_orig - solver.finiteDiff.y[i]) << std::endl;
    }*/
    return 0;

}
