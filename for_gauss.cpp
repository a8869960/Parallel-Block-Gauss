//
// Created by varsem on 15.11.23.
//
#include "functions.h"

#define eps 1e-15

int find_local_block_main(ARGS *arg, int step, int k, int l, double NORM)
{
    double min = 1.7976931348623158e+308, norm;
    int imax = step, jmax = step, count = 0;

    for(int i = step + arg->g; i < k; i += arg->p)
        for(int j = step; j < k; j++)
        {
            get_block(arg->a, arg->block, arg->indi[i], arg->indj[j], arg->n, arg->m, k, l);

            if (inverseMatrix(arg->block, arg->block_inv, arg->block_h, arg->m, arg->indi_m, arg->indj_m, NORM) == 0)
            {
                norm = matrixNorm(arg->block_inv, arg->m);

                if (norm < min)
                {
                    min = norm;

                    imax = i;
                    jmax = j;
                }
            } else
                count++;
        }

    if(count == (k - step) * (k - step))
        return -1;

    arg->local_imax = imax;
    arg->local_jmax = jmax;

    return 0;
}

int find_global_block_main(ARGS *arg, int k, int l, double NORM)
{
    int count = 0, imax = 0, jmax = 0;
    double norm = 0, min = 1.7976931348623158e+308;

    for(int i = -arg->g; i < arg->p - arg->g; i++)
    {
        if((arg + i)->status == io_status::no_matrix_main)
            count++;
        else
        {
            get_block(arg->a, arg->block, arg->indi[(arg + i)->local_imax], arg->indj[(arg + i)->local_jmax], arg->n, arg->m, k, l);

            inverseMatrix(arg->block, arg->block_inv, arg->block_h, arg->m, arg->indi_m, arg->indj_m, NORM);
            norm = matrixNorm(arg->block_inv, arg->m);

            if (norm < min)
            {
                min = norm;

                imax = (arg + i)->local_imax;
                jmax = (arg + i)->local_jmax;
            }
        }
    }
    if(count == arg->p)
        return -1;

    arg->global_imax = imax;
    arg->global_jmax = jmax;

    return 0;
}

void rearrange_elements(ARGS *arg, int step)
{
    int helper;

    helper = arg->indi[arg->global_imax];
    arg->indi[arg->global_imax] = arg->indi[step];
    arg->indi[step] = helper;

    helper = arg->indj[arg->global_jmax];
    arg->indj[arg->global_jmax] = arg->indj[step];
    arg->indj[step] = helper;
}

double matrixNorm(double *A, int n)
{
    double norm = 0, helper = 0;

    for(int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
            helper += abs(A[i * n + j]);

        if(helper > norm)
            norm = helper;

        helper = 0;
    }

    return norm;
}

void get_block(
        double* A,
        double* block,
        int i,
        int j,
        int n,
        int m,
        int k,
        int l)
{
    int block_m = (i == k ? l : m), block_l = (j == k ? l : m);

    int r, s;
    int a = i * n * m + j * m; //number of first element of the block
    memset(block, 0, sizeof(double) * m * m);

    for(r = 0; r < block_m; r++)
        for(s = 0; s < block_l; s++)
            block[r * block_l + s] = A[a + r * n + s];
}

void put_block(
        double* A,
        double* block,
        int i,
        int j,
        int n,
        int m,
        int k,
        int l)
{
    int block_m = (i == k ? l : m), block_l = (j == k ? l : m);

    int r, s;
    int a = i * n * m + j * m; //number of first element of the block

    for(r = 0; r < block_m; r++)
    {
        for(s = 0; s < block_l; s++)
        {
            A[a + r * n + s] = block[r * block_l + s];
        }
    }
}

void E(double* block, int m)
{
    memset(block, 0, sizeof(double) * m * m);

    for(int i = 0; i < m; i++)
        block[i * m + i] = 1;
}

void get_block_b( double *B, double *block, int i, int m, int k, int l)
{
    int block_m = (i == k ? l : m);

    memset(block, 0, sizeof(double) * m * m);

    int r;
    int b = i * m; //number of first element of the block

    for(r = 0; r < block_m; r++)
    {
        block[r] = B[b + r];
    }
}

void put_block_b( double *B, double *block, int i, int m, int k, int l)
{
    int block_m = (i == k ? l : m);

    int r;
    int b = i * m; //number of first element of the block

    for(r = 0; r < block_m; r++)
        B[b + r] = block[r];
}

void reduce_sum(int p, double* a, int n)
{
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static double *r = nullptr;
    int i;

    if(p <= 1)
        return;
    pthread_mutex_lock(&m);

    if(r == nullptr)
        r = a;
    else
        for(i = 0; i < n; i++) r[i] += a[i];

    t_in++;
    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else
        while(t_in < p)
            pthread_cond_wait(&c_in, &m);

    if(r != a)
        for(i = 0; i < n; i++) a[i] = r[i];

    t_out++;
    if(t_out >= p)
    {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    }
    else
        while(t_out < p)
            pthread_cond_wait(&c_out, &m);

    pthread_mutex_unlock(&m);
}

double get_full_time()
{
    struct timeval buf;
    gettimeofday(&buf, NULL);
    return buf.tv_sec + buf.tv_usec / 1.e6;
}

double get_CPU_time()
{
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1.e6;
}