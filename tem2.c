#include "stdio.h"

#include "temp.hpp"

int temp2_(int * nptr)
{
    temp(*nptr);
    printf("2!\n");

    return 0;
}