//
// Created by varsem on 15.11.23.
//
#include "functions.h"

void *process_gauss(void *arg_)
{
    ARGS *arg = (ARGS*)arg_;
    int g = arg->g, n = arg->n, m = arg->m, p = arg->p;

   int k, l, bl, step;
   double norm = matrixNorm(arg->A, n);

   double *a = arg->a, *b = arg->b, *block = arg->block, *block_h = arg->block_h, *block_inv = arg->block_inv;
   int *indi = arg->indi, *indj = arg->indj, *indi_m = arg->indi_m, *indj_m = arg->indj_m;

    k = n / m; //how many blocks m*m
    l = n - k * m; //how long last block
    bl = (l != 0) ? k + 1 : k; //number of all blocks

    double start_CPU, start_FULL, end_CPU, end_FULL;
    start_CPU = get_CPU_time();
    start_FULL = get_full_time();
    //Прямой ход
    for(step = 0; step < k; step++)
    {
        //Ищем главный элемент
        if(find_local_block_main(arg, step, k, l, norm) == -1)
        {
            arg->status = io_status::no_matrix_main;
        }
        reduce_sum(p);

        //Ищем "самый" главный
        if(find_global_block_main(arg, k, l, norm) == -1)
        {
            cout << "Can't find main block element on step " << step << endl;
            arg->status = io_status::no_matrix_main;
            end_CPU = get_CPU_time();
            end_FULL = get_full_time();
            arg->cpu_time = end_CPU - start_CPU;
            arg->full_time = end_FULL - start_FULL;
            return 0;
        }

        //Переставляем элементы
        if(g == 0)
            rearrange_elements(arg, step);
        reduce_sum(p);

        //Деление первой строчки
        get_block(a, block, indi[step], indj[step], n, m, k, l);
        inverseMatrix(block, block_inv, block_h, m, indi_m, indj_m, norm);

        for(int j = step + g + 1; j < bl; j += p)
        {
            int size_l = (j != k ? m : l);

            get_block(a, block, indi[step], indj[j], n, m, k, l);
            matrix_product(block_inv, block, block_h, m, m, size_l);
            put_block(a, block_h, indi[step], indj[j], n, m, k, l);
        }
        reduce_sum(p);

        if(g == 0)
        {
            //Делаем главный элемент единичкой
            E(block, m);
            put_block(a, block, indi[step], indj[step], n, m, k, l);

            //Делим присоединенный столбец
            get_block_b(arg->b, block, indi[step], m, k, l);
            matrix_product(block_inv, block, block_h, m, m, 1);
            put_block_b(arg->b, block_h, indi[step], m, k, l);
        }
        reduce_sum(p);

        //Зануление столбца
        for(int i = step + g + 1; i < bl; i += p)
        {
            get_block(a, block, indi[i], indj[step], n, m, k ,l); // size = size_m * m

            //Зануляем А
            for(int j = step + 1; j < bl; j++)
            {
                int size_m = (i != k ? m : l);
                int size_l = (j != k ? m : l);

                get_block(a, block_h, indi[step], indj[j], n, m, k, l); //size = m * size_l
                matrix_product(block, block_h, block_inv, size_m, m, size_l);
                get_block(a, block_h, indi[i], indj[j], n, m, k, l);
                matrixSubtraction(block_h, block_inv, block_h, size_m, size_l);
                put_block(a, block_h, indi[i], indj[j], n, m, k, l);
            }
            //Зануляем В
            get_block_b(b, block_h, indi[step], m, k, l);
            matrix_product(block, block_h, block_inv, m, m, 1);
            get_block_b(b, block_h, indi[i], m, k, l);
            matrixSubtraction(block_h, block_inv, block_h, 1, m);
            put_block_b(b, block_h, indi[i], m, k, l);
        }

        memset(block, 0, sizeof(double) * m * m);
        for(int i = step + g + 1; i < bl; i += p)
            put_block(a, block, indi[i], indj[step], n, m, k, l);
//        reduce_sum(p);
    }

    //Делим последний блок
    if(bl == k + 1 and g == 0)
    {
        get_block(a, block, indi[k], indj[k], n, m, k, l);
        if(inverseMatrix(block, block_inv, block_h, l, indi_m, indj_m, norm) == -1)
        {
            cout << "Can't find block max." << endl;
            arg->status = io_status::no_matrix_main;
        } else
        {
            E(block, l);
            put_block(a, block, indi[k], indj[k], n, m, k, l);

            get_block_b(b, block, indi[k], m, k, l);
            matrix_product(block_inv, block, block_h, l, l, 1);
            put_block_b(b, block_h, indi[k], m, k, l);
        }
    }
    reduce_sum(p);

    if((arg - arg->g)->status == io_status::no_matrix_main)
    {
        end_CPU = get_CPU_time();
        end_FULL = get_full_time();
        arg->cpu_time = end_CPU - start_CPU;
        arg->full_time = end_FULL - start_FULL;
        return 0;
    }

    //Обратный ход
    if(bl == k + 1)
    {
        get_block_b(b, block, indi[k], m, k, l); //size = 1 * l or l * 1

        for(int i = k - g - 1; i >= 0; i -= p)
        {
            get_block(a, block_h, indi[i], indj[k], n, m, k, l); //size = m * l
            matrix_product(block_h, block, block_inv, m, l, 1); //size = m * 1 or 1 * m
            get_block_b(b, block_h, indi[i], m, k, l); //size = 1 * m or m * 1
            matrixSubtraction(block_h, block_inv, block_h, 1, m);
            put_block_b(b, block_h, indi[i], m, k, l);
        }
    }
    reduce_sum(p);

    for(step = k - 1; step >= 0; step--)
    {
        get_block_b(b, block, indi[step], m, k, l); //size = 1 * m or m * 1

        for(int i = step - g - 1; i >= 0; i -= p)
        {
            get_block(a, block_h, indi[i], indj[step], n, m, k, l); //size = m * m
            matrix_product(block_h, block, block_inv, m, m, 1); //size = m * 1 or 1 * m
            get_block_b(b, block_h, indi[i], m, k, l); //size = 1 * or m * 1
            matrixSubtraction(block_h, block_inv, block_h, 1, m); //
            put_block_b(b, block_h, indi[i], m, k, l);
        }
    }
    reduce_sum(p);

    end_CPU = get_CPU_time();
    end_FULL = get_full_time();

    if(g == 0)
    {
        for (int i = 0; i < k; i++)
            for (int j = 0; j < m; j++)
                arg->x[indj[i] * m + j] = b[indi[i] * m + j];
        for (int j = 0; j < l; j++)
            arg->x[m * k + j] = b[indi[k] * m + j];
    }

    arg->status = io_status::success;
    arg->cpu_time = end_CPU - start_CPU;
    arg->full_time = end_FULL - start_FULL;

    return nullptr;
}