//
// Created by Лаэтин  Андрей on 06.01.2023.
//

#include "Solver.h"

Solver::Solver(int n, double a, double b){
    finiteDiff = FDM(n, a, b);
    galerkinMethod = GalerkinMethod(a, b ,n);
}

void Solver::SolveFDM(int type){
    finiteDiff.solve(type);
}
void Solver::SolveGal(){
    this->galerkinMethod.solve();
}
