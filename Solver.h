//
// Created by Лаэтин  Андрей on 06.01.2023.
//

#ifndef MDA_SOLVER_H
#define MDA_SOLVER_H

#include "FDM.h"
#include "GalerkinMethod.h"

class Solver{
public:
    Solver(int n, double a, double b);
    FDM finiteDiff;
    GalerkinMethod galerkinMethod;

    void SolveGal();
    void SolveFDM(int type);
};
#endif //MDA_SOLVER_H
