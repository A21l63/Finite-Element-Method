#include <iostream>
#include <cmath>
#include <fstream>

//Fk = h*h*fk
double f_K(double x) {
    return exp(x);//(- M_PI * M_PI) * x * cos(M_PI * x); //- 2 * M_PI * M_PI * sin(M_PI * x);
}

double p_k(double x) {
    return (exp(x) + 1);
}

double F_K(double x, double h) {
    return h * h * f_K(x);
}

double r_K(double x) {
    return -exp(x);
}

double B_K(double x, double h) {
    return 2 * p_k(x) - r_K(x) * h * h;
}


double q_K(double x) {
    return -1;
}

double A_K(double x, double h) {
    return p_k(x) + h * q_K(x) / 2;
}

// Ck = 1 - q_k(h/2)

double C_K(double x, double h) {
    return p_k(x) - h * q_K(x) / 2;
}

    //lambda_K
void lambda_K(double* lambda, double* X, double a_0, double a_1, double b_0, double b_1, int n) {
    double h = X[1] - X[0];
    lambda[0] = (a_1 * (B_K(X[1], h)/ A_K(X[1], h) - 4) / (2 * h)) / (a_0 + a_1 * (C_K(X[1], h) / A_K(X[1], h) - 3) / (2 * h));
    //lambda[0] = (-a_1 / h) / (a_0 - a_1 / h);
    // 1 ~ n-2
    for (int i = 1; i < n; i++) {

        lambda[i] = A_K(X[i], h) / (B_K(X[i], h) - C_K(X[i], h) * lambda[i - 1]) ;
    }
    //  n - 1
    lambda[n] = - ((b_1 / (2 * h)) * (B_K(X[n - 1], h) / C_K(X[n - 1], h) - 4)) / (b_0 + b_1 * (3 - A_K(X[n - 1], h) / C_K(X[n - 1], h)) / (2 * h));
    //lambda[n] = (b_1 / h) / (b_0 + b_1 / h);
}
    //gamma_K

void gamma_K(double* gamma, double* X, double* lambda, double a_0, double a_1, double b_0, double b_1, double A, double B, int n) {
    double h = X[1] - X[0];

    gamma[0] = (A + a_1 * F_K(X[1], h) / (A_K(X[1], h) * 2 * h)) / (a_0 + a_1 * (C_K(X[1], h) / A_K(X[1], h) - 3) / (2 * h));
    //gamma[0] = A / (a_0 - a_1 / h);
    // 1 ~ n - 2
    for (int i = 1; i < n; i++) {
        gamma[i] = - (F_K(X[i], h) - C_K(X[i], h) * gamma[i - 1])/(B_K(X[i], h) - C_K(X[i], h) * lambda[i - 1]);
    }

    gamma[n] = (B - b_1 * F_K(X[n - 1], h) / (C_K(X[n - 1], h) * 2 * h)) / (b_0 + b_1 * (3 - A_K(X[n - 1], h) / C_K(X[n - 1], h)) / (2 * h));
    //gamma[n] = B / (b_0 + b_1 / h);
}

double y_n(double lambda_n, double gamma_n, double lambda_n1, double gamma_n1) {
    return (lambda_n * gamma_n1 + gamma_n) / (1 - lambda_n * lambda_n1);
}

void grid(double* X, int n, double a, double b) {
    double h = (b - a) / n;
    for (int i = 1; i < n; i++) {
        X[i] = a + i * h;
    }
    X[0] = a;
    X[n] = b;
}


int main() {
        int n = 10000;
        double a = 0;
        double b = 1;
        double A = -3;
        double B = -1;
        double a_0 = 2;
        double a_1 = -3;
        double b_0 = 1;
        double b_1 = -1;


        double X[n + 1];
        grid(X, n, a, b);
        double y[n + 1];
        double lambda[n + 1];
        double gamma[n + 1];


        lambda_K(lambda, X, a_0, a_1, b_0, b_1, n);
        gamma_K(gamma, X, lambda, a_0, a_1, b_0, b_1, A, B, n);

        y[n] = (lambda[n] * gamma[n - 1] + gamma[n]) / (1 - lambda[n] * lambda[n - 1]);
        for (int i = n - 1; i >= 0; i--) {
            y[i] = lambda[i] * y[i + 1] + gamma[i];
        }

        double max = 1000;
        double y_orig;std::cout << "x |  y[i] | sin(PI * x) | * \n";
        for (int i = 0; i <= n; i++) {
            y_orig = exp(X[i]) - 1; //X[i] * cos(M_PI * X[i]); //sin(M_PI * X[i]);
            std::cout << X[i] << " | " << y[i] << " | " << y_orig << " | " << std::abs(y_orig - y[i]) << std::endl;
        }
        std::ofstream out;
        out.open("/Users/laetinandrej/CLionProjects/mda/output100.txt");
         if (out.is_open()) {
             for (int i = 0; i <= n; i++) {
                 y_orig = exp(X[i]) - 1; //sin(M_PI * X[i]);
                 out << X[i] << " , " << y[i] << " , " << y_orig << " , " << std::abs(y_orig - y[i]) << std::endl;
             }
         }
    return 0;

}