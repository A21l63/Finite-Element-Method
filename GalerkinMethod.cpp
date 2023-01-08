//
// Created by Лаэтин Андрей on 08.01.2023.
//

#include <iostream>
#include "GalerkinMethod.h"

GalerkinMethod::GalerkinMethod(double a, double b, int num) {
    numberOfNodes = num;
    numberOfFunc = numberOfNodes - 2;
    grid = Grid(numberOfNodes - 1, a, b);
    h = (b - a) / numberOfNodes;
    integralStep = 2 * h / intAcc;

    result.resize(numberOfFunc);
    A.resize(numberOfFunc, -1/h);
    B.resize(numberOfFunc,  2/h);
    C.resize(numberOfFunc, -1/h);
    A[0] = 0;
    C[numberOfFunc] = 0;
    D.resize(numberOfFunc);

    for (int i = 0; i < numberOfFunc; i++){
        D[i] = trapInt(grid.X[i], grid.X[i + 2]);
    }
}

double GalerkinMethod::funct(double a, double b, double x){
    if (x <= (a + h)){
        return (x - a) / h;
    }else{
        return (b - x) / h;
    }

}

double GalerkinMethod::trapInt(double a, double b) {
    double simpson_integral = 0;
    for(int step = 0; step < intAcc; step++) {
        double x1 = a + step * integralStep;
        double x2 = a + (step+1) * integralStep;

        simpson_integral += h / 6.0 * (equation.f_K(x1) * funct(a, b, x1) +
                4.0 * equation.f_K(0.5 * (x1 + x2)) * funct(a, b, 0.5 * (x1 + x2)) +
                equation.f_K(x2) * funct(a, b, x2));
    }
    return simpson_integral;
}

void GalerkinMethod::triDiag(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d, int n){
    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }

    result[0] = 0;
    result [numberOfNodes - 1] = 0;
    for (int i = 0; i < d.size(); ++i) {
        result[i + 1] = d[i];
    }
}


GalerkinMethod::GalerkinMethod(){}

double GalerkinMethod::L2norm(){
    double simpson_integral = 0;
    for (int i = 0; i < numberOfNodes; i += 2) {
        simpson_integral += h / 6.0 * (((grid.X[i]) * (1 - grid.X[i]) / 2 - result[i]) +
                4.0 * ((grid.X[i + 1] *(1 - grid.X[i + 1]))/2 - result[i + 1])+
                (grid.X[i + 2]) * (1 - grid.X[i + 2]) / 2 - result[i + 2]);
        }
    return simpson_integral;
}

void GalerkinMethod::solve() {
    double y_orig;
    triDiag(this->A,this->B,this->C, this->D, numberOfFunc);
    double norm = 0;
    double abbs = 0;
    for (int i = 0; i < numberOfNodes; i++) {
        y_orig = (grid.X[i]) * (1 - grid.X[i])/2;
        abbs = abs(result[i] - y_orig);
        //std::cout << "x = " << grid.X[i] << " Result = " << result[i] << " Correct = " << y_orig << " ||u - y || " << abs(result[i] - y_orig) << std::endl;
        if (norm <= abbs) {
            norm = abbs;
        }
    }
    std::cout << L2norm();
}
