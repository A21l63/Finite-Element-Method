//
// Created by Лаэтин  Андрей on 26.11.2022.
//

#ifndef MDA_EQUATION_H
#define MDA_EQUATION_H
#include <cmath>


class Equation {
public:
    double A = 0;
    double B = 0;
    double a_0 = 1;
    double a_1 =0;
    double b_0 = 1;
    double b_1 = 0;

    double f_K(double x) {
        return 2-3 * x;//(- M_PI * M_PI) * x * cos(M_PI * x); //- 2 * M_PI * M_PI * sin(M_PI * x);
    }

    double p_k(double x) {
        return 1;
    }

    double F_K(double x, double h) {
        return h * h * f_K(x);
    }

    double r_K(double x) {
        return 0;
    }

    double B_K(double x, double h) {
        return 2 * p_k(x) - r_K(x) * h * h;
    }


    double q_K(double x) {
        return 0;
    }

    double A_K(double x, double h) {
        return p_k(x) + h * q_K(x) / 2;
    }

// Ck = 1 - q_k(h/2)

    double C_K(double x, double h) {
        return p_k(x) - h * q_K(x) / 2;
    }
};


#endif //MDA_EQUATION_H
