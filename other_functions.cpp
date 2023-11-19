//
// Created by varsem on 15.11.23.
//
#include <stdlib.h>
#include <errno.h>
#include <limits.h>

#include "functions.h"

int toInt(const char* str, int* ptr)
{
    long L;
    char* e;

    errno = 0;
    L = strtol(str, &e, 10);

    if (!errno && *e == '\0')
        if (INT_MIN <= L && L <= INT_MAX)
        {
            *ptr = (int)L;
            return 0;
        }
        else
            return -1;
    else
        return -1;
}
