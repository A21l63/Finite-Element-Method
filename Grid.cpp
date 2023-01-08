//
// Created by Лаэтин  Андрей on 26.11.2022.
//

#include "Grid.h"

Grid::Grid(int n, double a, double b) {
    X = new double(n+1);
    num = n;

    h = (b - a) / n;

    for (int i = 1; i < n; i++) {
        X[i] = a + i * h;
    }
    X[0] = a;
    X[n] = b;
}

Grid::Grid() {}

