//
// Created by varsem on 15.11.23.
//
#include <iostream>
#include <cstring>
#include <pthread.h>
#include <ctime>

#include <cstdio>
#include <cstdlib>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

enum class io_status
{
    success,
    none,
    no_matrix_main
};

class ARGS
{
public:
    int n = 0;
    int m = 0;
    double *A = nullptr;
    double *B = nullptr;
    double *x = nullptr;
    int *indi_m = nullptr;
    int *indj_m = nullptr;
    int *indi = nullptr;
    int *indj = nullptr;
    double *a = nullptr;
    double *b = nullptr;
    double *block = nullptr;
    double *block_inv = nullptr;
    double *block_h = nullptr;

    int g = 0;
    int p = 0;

    io_status status = io_status::none;

    int local_imax = 0;
    int local_jmax = 0;
    int global_imax = 0;
    int global_jmax = 0;

    double cpu_time = 0;
    double full_time = 0;
    
    ~ARGS()
    {
        delete[] indi_m;
        delete[] indj_m;
        delete[] block;
        delete[] block_h;
        delete[] block_inv;
    }

    void print()
    {
        cout << "m      | " << g << endl;
        cout << "cpu    | " << cpu_time << endl;
        cout << "full   | " << full_time << endl;
        cout << "status | ";
        switch (status)
        {
            case(io_status::success):
                cout << "success" << endl;
                break;

            case(io_status::no_matrix_main):
                cout << "no_matrix_main" << endl;
                break;

            case(io_status::none):
                cout << "none" << endl;
                break;
        }
    }
};

//other_functions
int toInt(const char* str, int* ptr);

//matrix_input
int fileMatrixInput(double* A, char* filename, int n);
void f(double *A, int s, int n);
void init_B(double *B, double *A, int n);

//matrix_output
void matrixOutput(double *matrix, int l, int n, int r);

//results
int calc_r1(double* A, double* x, double* B, int n, double* helper, double *r1);
void calc_r2(double *x, int n, double *r2);

//matrix_operations
int matrixSubtraction(double* A1, double *A2, double *C, int n, int m);
void matrix_product(double *A, double* B, double* C, int n, int s, int m);

//process_gauss
void *process_gauss(void *arg);

//inverse_matrix
int inverseMatrix(double *a, double *A, double *B, int n, int *indi_m, int *indj_m, double norm); //a^(-1) = A

//for_gauss
int find_local_block_main(ARGS *arg, int step, int k, int l, double NORM);
int find_global_block_main(ARGS *arg, int k, int l, double NORM);
void rearrange_elements(ARGS *arg, int step);

double matrixNorm(double *A, int n);

void get_block(
        double* A,
        double* block,
        int i,
        int j,
        int n,
        int m,
        int k,
        int l);
void put_block(
        double* A,
        double* block,
        int i,
        int j,
        int n,
        int m,
        int k,
        int l);

void E(double* block, int m);

void get_block_b( double *B, double *block, int i, int m, int k, int l);
void put_block_b( double *B, double *block, int i, int m, int k, int l);

void reduce_sum(int p, double* a = nullptr, int n = 0);

double get_full_time();
double get_CPU_time();

