//
// Created by varsem on 15.11.23.
//
#include "functions.h"
#include "gauss.h"

#include <cstring>

int gauss_func(int n,
               int m,
               double *A,
               double *B,
               double *x,
               int *indi_m,
               int *indj_m,
               int *indi,
               int *indj,
               double* a,
               double *b,
               double *block,
               double *block_inv,
               double *block_h)
{
    int k, l, bl, i, j, step;
    double norm = matrixNorm(A, n);

    k = n / m; //how many blocks m*m
    l = n - k * m; //how long last block
    bl = (l != 0) ? k + 1 : k; //number of all blocks

    for(i = 0; i < bl; i++)
    {
        indi[i] = i;
        indj[i] = i;
    }

    memcpy(a, A, sizeof(double) * n * n);
    memcpy(b, B, sizeof(double) * n);

    //Прямой ход
    for(step = 0; step < k; step++)
    {
        if(matrixMax(a, step, n, m, k, l, indi, indj, indi_m, indj_m, block, block_inv, block_h, norm) == -1)
        {
            cout << "Can't find block max on step " << step << endl;
            return -1;
        }

        get_block(a, block, indi[step], indj[step], n, m, k, l);

//        cout << "BLOCK II DO " << indi[step] << indj[step] << endl;
//        matrixOutput(block, m, m, 12);

        inverseMatrix(block, block_inv, block_h, m, indi_m, indj_m, norm);

//        cout << "BLOCK II" << endl;
//        matrixOutput(block_inv, m, m, 12);
//        cout << endl;
//        cout << "MATRIX PR" << endl;
//        matrix_product(block_inv, block, block_h, m, m, m);
//        matrixOutput(block_h, m, m, 12);
//
        //Деление первой строчки
        E(block, m);
        put_block(a, block, indi[step], indj[step], n, m, k, l);
        for(j = step + 1; j < bl; j++)
        {
            int size_l = (j != k ? m : l);

            get_block(a, block, indi[step], indj[j], n, m, k, l);
//            cout << "BLOCK" << endl;
//            matrixOutput(block, m, size_l, 12);
            matrix_product(block_inv, block, block_h, m, m, size_l);
//            cout << "BLOCK H" << endl;
//            matrixOutput(block_h, m, size_l, 12);
            put_block(a, block_h, indi[step], indj[j], n, m, k, l);
        }
        get_block_b(b, block, indi[step], m, k, l);
//        cout << "B !!!" << endl;
//        matrixOutput(block, m, 1, 12);
        matrix_product(block_inv, block, block_h, m, m, 1);
//        cout << "BH !!!" << endl;
//        matrixOutput(block_h, m, 1, 12);
        put_block_b(b, block_h, indi[step], m, k, l);

//        cout << "A /" << step + 1<< endl;
//        matrixOutput(a, n, n, 12);
//        cout << "B /" << endl;
//        matrixOutput(b, 1, n, 12);
//        cout << endl;

        //Зануление столбца
        for(i = step + 1; i < bl; i++)
        {
            get_block(a, block, indi[i], indj[step], n, m, k ,l); // size = size_m * m
//            cout << "BLOCK" << endl;
//            matrixOutput(block, m, m, 3);
            for(j = step + 1; j < bl; j++)
            {
                int size_m = (i != k ? m : l);
                int size_l = (j != k ? m : l);

                get_block(a, block_h, indi[step], indj[j], n, m, k, l); //size = m * size_l

//                cout << "BLOCK" << endl;
//                matrixOutput(block, size_m, m, 12);
//                cout << "BLOCK_H" << endl;
//                matrixOutput(block_h, m, size_l, 12);

                matrix_product(block, block_h, block_inv, size_m, m, size_l);

//                cout << "BLOCK_INV" << endl;
//                matrixOutput(block_inv, size_m, size_l, 12);

                get_block(a, block_h, indi[i], indj[j], n, m, k, l);
//                cout << "BLOCK_H" << endl;
//                matrixOutput(block_h, size_m, size_l, 3);
                matrixSubtraction(block_h, block_inv, block_h, size_m, size_l);
//                cout << "BLOCK_H" << endl;
//                matrixOutput(block_h, size_m, size_l, 3);
                put_block(a, block_h, indi[i], indj[j], n, m, k, l);
            }
            get_block_b(b, block_h, indi[step], m, k, l);
//            cout << "BLOCK_H" << endl;
//            matrixOutput(block_h, m, m, 3);
            matrix_product(block, block_h, block_inv, m, m, 1);

            get_block_b(b, block_h, indi[i], m, k, l);
            matrixSubtraction(block_h, block_inv, block_h, 1, m);
            put_block_b(b, block_h, indi[i], m, k, l);
        }
        memset(block, 0, sizeof(double) * m * m);
        for(i = step + 1; i < bl; i++)
            put_block(a, block, indi[i], indj[step], n, m, k, l);

//        cout << "A -" << endl;
//        matrixOutput(a, n, n, 12);
//        cout << "B -" << endl;
//        matrixOutput(b, 1, n, 12);
//        cout << endl;
    }

//    cout << "PR HOD A1" << endl;
//    matrixOutput(a, n, n, 12);
//    cout << "PR HOD B1" << endl;
//    matrixOutput(b, 1, n, 12);
//    cout << endl;

    if(bl == k + 1)
    {
        get_block(a, block, indi[k], indj[k], n, m, k, l);
//        cout << "BLOCK LL" << endl;
//        matrixOutput(block, l, l, 12);
        if(inverseMatrix(block, block_inv, block_h, l, indi_m, indj_m, norm) == -1)
        {
            cout << "Can't find block max." << endl;
            return -1;
        }

//        cout << "BLOCK II" << endl;
//        matrixOutput(block_inv, l, l, 7);

        E(block, l);
//        cout << "BLOCK LL" << endl;
//        matrixOutput(block, l, l, 12);
        put_block(a, block, indi[k], indj[k], n, m, k, l);
//        cout << "PR HOD A" << endl;
//        matrixOutput(a, n, n, 5);

        get_block_b(b, block, indi[k], m, k, l);
        matrix_product(block_inv, block, block_h, l, l, 1);
        put_block_b(b, block_h, indi[k], m, k, l);
    }

//   cout << "PR HOD A" << endl;
//   matrixOutput(a, n, n, 12);
//   cout << "PR HOD B" << endl;
//   matrixOutput(b, 1, n, 12);
//   cout << endl;

    //Обратный ход
    if(bl == k + 1)
    {
        get_block_b(b, block, indi[k], m, k, l); //size = 1 * l or l * 1
//        cout << "BLOCK" << endl;
//        matrixOutput(block, m, m, 3);

        for(i = k - 1; i >= 0; i--)
        {
            get_block(a, block_h, indi[i], indj[k], n, m, k, l); //size = m * l
//            cout << "BLOCK_H" << endl;
//            matrixOutput(block_h, m, m, 3);

            matrix_product(block_h, block, block_inv, m, l, 1); //size = m * 1 or 1 * m
//            cout << "BLOCK_INV" << endl;
//            matrixOutput(block_inv, m, m, 3);

            get_block_b(b, block_h, indi[i], m, k, l); //size = 1 * m or m * 1
//            cout << "BLOCK_H" << endl;
//            matrixOutput(block_h, 1, m, m);

            matrixSubtraction(block_h, block_inv, block_h, 1, m);
//            cout << "BLOCK_H" << endl;
//            matrixOutput(block_h, 1, m, m);

            put_block_b(b, block_h, indi[i], m, k, l);
        }
    }
    for(step = k - 1; step >= 0; step--)
    {
        get_block_b(b, block, indi[step], m, k, l); //size = 1 * m or m * 1

        for(i = step - 1; i >= 0; i--)
        {
            get_block(a, block_h, indi[i], indj[step], n, m, k, l); //size = m * m
            matrix_product(block_h, block, block_inv, m, m, 1); //size = m * 1 or 1 * m
            get_block_b(b, block_h, indi[i], m, k, l); //size = 1 * or m * 1
            matrixSubtraction(block_h, block_inv, block_h, 1, m); //
            put_block_b(b, block_h, indi[i], m, k, l);
        }
    }

//    cout << "OB HOD A" << endl;
//    matrixOutput(a, n, n, 12);
//    cout << "OB HOD B" << endl;
//    matrixOutput(b, 1, n, 12);
//    cout << endl;

    for(i = 0; i < k; i++)
        for(j = 0; j < m; j++)
            x[indj[i] * m + j] = b[indi[i] * m + j];
    for(j = 0; j < l; j++)
        x[m * k + j] = b[indi[k] * m + j];
//    memcpy(x, b, sizeof(double) * n);

    return 0;
}
