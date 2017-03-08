#include "stdinc.h"

double xrandom(double xl, double xh)
{

    return (xl + (xh - xl) * ((double) rand()) / (double)RAND_MAX);
}


