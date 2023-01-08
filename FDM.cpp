//
// Created by Лаэтин  Андрей on 06.01.2023.
//

#include "FDM.h"

FDM::FDM(int n, double a, double b) {
    grid = Grid(n, a, b);
    y = new double[n+1];
    lambda = new double[n+1];
    gamma = new double[n+1];
}

FDM::FDM(){}

void FDM::lambda_K_1() {
    lambda[0] = (-equation.a_1 / grid.h) / (equation.a_0 - equation.a_1 / grid.h);
    for (int i = 1; i < grid.num; i++) {
        lambda[i] = equation.A_K(grid.X[i], grid.h) / (equation.B_K(grid.X[i], grid.h) - equation.C_K(grid.X[i], grid.h) * lambda[i - 1]) ;
    }
    lambda[grid.num] = (equation.b_1 / grid.h) / (equation.b_0 + equation.b_1 / grid.h);
}

void FDM::gamma_K_1() {
    gamma[0] = equation.A / (equation.a_0 - equation.a_1 / grid.h);
    for (int i = 1; i < grid.num; i++) {
        gamma[i] = -(equation.F_K(grid.X[i], grid.h) - equation.C_K(grid.X[i], grid.h) * gamma[i - 1]) /
                   (equation.B_K(grid.X[i], grid.h) - equation.C_K(grid.X[i], grid.h) * lambda[i - 1]);
    }
    gamma[grid.num] = equation.B / (equation.b_0 + equation.b_1 / grid.h);
}

void FDM::lambda_K_2() {
    lambda[0] = (equation.a_1 * (equation.B_K(grid.X[1], grid.h)/ equation.A_K(grid.X[1], grid.h) - 4) / (2 * grid.h)) / (equation.a_0 + equation.a_1 * (equation.C_K(grid.X[1], grid.h) / equation.A_K(grid.X[1], grid.h) - 3) / (2 * grid.h));

    for (int i = 1; i < grid.num; i++) {
        lambda[i] = equation.A_K(grid.X[i], grid.h) / (equation.B_K(grid.X[i], grid.h) - equation.C_K(grid.X[i], grid.h) * lambda[i - 1]) ;
    }
    lambda[grid.num] = - ((equation.b_1 / (2 * grid.h)) * (equation.B_K(grid.X[grid.num - 1], grid.h) / equation.C_K(grid.X[grid.num - 1], grid.h) - 4)) / (equation.b_0 + equation.b_1 * (3 - equation.A_K(grid.X[grid.num - 1], grid.h) / equation.C_K(grid.X[grid.num - 1], grid.h)) / (2 * grid.h));
}

void FDM::gamma_K_2() {
    gamma[0] = (equation.A + equation.a_1 * equation.F_K(grid.X[1], grid.h) / (equation.A_K(grid.X[1], grid.h) * 2 * grid.h)) / (equation.a_0 + equation.a_1 * (equation.C_K(grid.X[1], grid.h) / equation.A_K(grid.X[1], grid.h) - 3) / (2 * grid.h));

    for (int i = 1; i < grid.num; i++) {
        gamma[i] = -(equation.F_K(grid.X[i], grid.h) - equation.C_K(grid.X[i], grid.h) * gamma[i - 1]) /
                   (equation.B_K(grid.X[i], grid.h) - equation.C_K(grid.X[i], grid.h) * lambda[i - 1]);
    }
    gamma[grid.num] = (equation.B - equation.b_1 * equation.F_K(grid.X[grid.num - 1], grid.h) / (equation.C_K(grid.X[grid.num - 1], grid.h) * 2 * grid.h)) / (equation.b_0 + equation.b_1 * (3 - equation.A_K(grid.X[grid.num - 1], grid.h) / equation.C_K(grid.X[grid.num - 1], grid.h)) / (2 * grid.h));
}

void FDM::solve(int type) {
    if (!type){
        lambda_K_1();
        gamma_K_1();
    } else
    {
        lambda_K_2();
        gamma_K_2();
    }

    y[grid.num] = (lambda[grid.num] * gamma[grid.num - 1] + gamma[grid.num]) / (1 - lambda[grid.num] * lambda[grid.num - 1]);
    for (int i = grid.num - 1; i >= 0; i--) {
        y[i] = lambda[i] * y[i + 1] + gamma[i];
    }

}



