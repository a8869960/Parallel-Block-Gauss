//
// Created by varsem on 15.11.23.
//
#include "functions.h"

#include <cstring>

double norm1(double *x, int n)
{
    double norma = 0;

    for(int i = 0; i < n; i++)
        norma += abs(x[i]);

    return norma;
}

int calc_r1(double* A, double* x, double* B, int n, double* helper, double *r1)
{
    memset(helper, 0, n);

    matrix_product(A, x, helper, n, n, 1);

    if(matrixSubtraction(helper,B, helper, n, 1) == -1)
        return -1;

    double help = norm1(B, n);

    if(abs(help) < 1e-16)
    {
        cout << "Division by 0." << endl;
        return 0;
    }

    *r1 = norm1(helper, n) / help;

    return 0;
}

void calc_r2(double *x, int n, double *r2)
{
    double s = 0;
    for(int i = 0; i < n; i++)
    {
        s += abs(x[i] - ((i + 1) % 2));
    }

    *r2 = s;
}