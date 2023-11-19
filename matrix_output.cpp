//
// Created by varsem on 15.11.23.
//
#include <cstdio>

#include "functions.h"

int min(int r, int l);

void matrixOutput(double *matrix, int l, int n, int r)
{
    int minN = min(r, l), minM = min(r, n);

    for(int i = 0; i < minN; i++)
    {
        for(int j = 0; j < minM; j++)
        {
            printf(" %10.3e", matrix[i * n + j]);
        }
        printf("\n");
    }
}

int min(int r, int l)
{
    return (r < l) ? r : l;
}