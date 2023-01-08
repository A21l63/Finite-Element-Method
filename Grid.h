//
// Created by Лаэтин  Андрей on 26.11.2022.
//

#ifndef MDA_GRID_H
#define MDA_GRID_H


class Grid {
public:
    int num;
    double h;

    double* X;

    Grid(int n, double a, double b);

    Grid();

};


#endif //MDA_GRID_H
