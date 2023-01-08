//
// Created by Лаэтин  Андрей on 06.01.2023.
//

#ifndef MDA_FDM_H
#define MDA_FDM_H
#include "Grid.h"
#include "Equation.h"

class FDM{
public:
    Grid grid;
    Equation equation;

    double* y;
    double* lambda;
    double* gamma;

    FDM(int n, double a, double b);

    FDM();


    void lambda_K_1();
    void gamma_K_1();
    void lambda_K_2();
    void gamma_K_2();

    void solve(int type);

};
#endif //MDA_FDM_H
