#include <math.h>

double f(double x, double y)
{
    return 2*((1+x)*sin(x+y)-cos(x+y));
}

double boundaries(double x, double y)
{
    return (1+x)*sin(x+y);
}
