//
// Created by Лаэтин  Андрей on 08.01.2023.
//

#ifndef MDA_GALERKINMETHOD_H
#define MDA_GALERKINMETHOD_H

#include <memory>
#include <vector>
#include "Equation.h"
#include "Grid.h"

class GalerkinMethod {
private:
    int intAcc = 2;
    double integralStep;
    int numberOfFunc;
    int numberOfNodes;
    double h;

    std::vector<double> result;
    std::vector<double> A;
    std::vector<double> B;
    std::vector<double> C;
    std::vector<double> D;
    double L2norm();

    Equation equation;
    Grid grid;
    double trapInt(double a, double b);
    void triDiag(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d, int n);
    double funct(double a, double b, double x);

public:
    GalerkinMethod(double a, double b, int num);
    GalerkinMethod();
    void solve();
};


#endif //MDA_GALERKINMETHOD_H
