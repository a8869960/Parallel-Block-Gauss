//
// Created by varsem on 15.11.23.
//
#include "functions.h"

int main(int ac, char* av[])
{
    try
    {
        int task = 11;

        int n, m, p, r, s;
        char* filename = nullptr;
        double *A = nullptr;

        if(ac != 6 and ac != 7)
        {
            cout << "Wrong parameters." << endl;
            return -1;
        }

        if(toInt(av[1], &n) == -1)
            return -1;

        if(toInt(av[2], &m) == -1)
            return -1;

        if(toInt(av[3], &p) == -1)
            return -1;

        if(toInt(av[4], &r) == -1)
            return -1;

        if(toInt(av[5], &s) == -1)
            return -1;

        if(s == 0 and ac == 7)
        {
            filename = av[6];
        } else if((s == 0 and ac == 6) or (s != 0 and ac == 7))
        {
            cout << "Wrong parameter s and filename." << endl;
            return -1;
        }

        //Проверка аргументов командной строки
        if(n < 1 or m < 1 or n < m)
        {
            cout << "Wrong n and m parameters." << endl;
            return -1;
        }
        if(s < 0 or s > 4)
        {
            cout << "Wrong parameter s." << endl;
            return -1;
        }
        if(r < 1)
        {
            cout << "Wrong parameter r." << endl;
            return -1;
        }

        //Инициализация матрицы A
        A = new double[n*n];

        if(ac == 6)
        {
            f(A, s, n);
        } else
        {
            if(fileMatrixInput(A, filename, n) == -1)
            {
                delete[] A;
                return -1;
            }
        }

        cout << "Matrix A:" << endl;
        matrixOutput(A, n, n, r);

        double *B = new double[n];
        init_B(B, A, n);

        cout << "Matrix B:" << endl;
        matrixOutput(B, 1, n, r);

        double *x = new double[n];
        double *helper = nullptr;

        //Метод Гаусса
        int k = n / m; //how many blocks m*m
        int l = n - k * m; //how long last block
        int bl = (l != 0) ? k + 1 : k; //number of all blocks

        bool flag = true;

        double *a = new double[n * n], *b = new double[n];
        int *indi = new int[bl], *indj = new int[bl];

        ARGS *args = new ARGS[p];
        pthread_t *threads = new pthread_t[p];

        //Заполняем аргументы
        for(int i = 0; i < bl; i++)
        {
            indi[i] = i;
            indj[i] = i;
        }

        memcpy(a, A, sizeof(double) * n * n);
        memcpy(b, B, sizeof(double) * n);

        for(int i = 0; i < p; i++)
        {
            args[i].n = n;
            args[i].m = m;
            args[i].A = A;
            args[i].B = B;
            args[i].x = x;
            args[i].indi_m = new int[m];
            args[i].indj_m = new int[m];
            args[i].indi = indi;
            args[i].indj = indj;
            args[i].a = a;
            args[i].b = b;
            args[i].block = new double[m * m];
            args[i].block_inv = new double[m * m];
            args[i].block_h = new double[m * m];

            args[i].g = i;
            args[i].p = p;
        }

        //Запускаем потоки
        clock_t start_time =  clock();
        for(int i = 0; i < p; i++)
        {
            if(pthread_create(threads + i, 0, process_gauss, args + i))
            {
                cout << "Cannot create thread " << i << endl;
                delete[] threads;
                delete[] args;
                return -4;
            }
        }

        //Ждем потоки
        for(int i = 0; i < p; i++)
        {
            if(pthread_join(threads[i], 0))
                cout << "Cannot wait thread " << i << endl;
        }

        int count = 0;
        for(int i = 0; i < p; i++)
            if(args[i].status != io_status::success)
                count++;

        if(count != 0)
        {
            cout << "Can't be used this method." << endl;
            flag = false;
        } else
        {
            cout << "Solution:" << endl;
            matrixOutput(x, 1, n, r);
        }

        //Запись результата
        double t1 = args[0].full_time;
        start_time = clock();
        double r1 = 0, r2 = 0;
        if(flag)
        {
            helper = new double[n];
            if (calc_r1(A, x, B, n, helper, &r1) == -1)
                return -1;

            calc_r2(x, n, &r2);
        }
        double t2 = (clock() - start_time) / CLOCKS_PER_SEC;

        printf (
                "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
                av[0], task, r1, r2, t1, t2, s, n, m, p);

        delete[] A;
        delete[] B;
        delete[] x;
        delete[] helper;
        delete[] a;
        delete[] b;
        delete[] threads;
        delete[] args;
        delete[] indi;
        delete[] indj;

        return 0;
    } catch (const bad_alloc& e)
    {
        cout << "Bad alloc" << endl;
        return -1;
    }
}