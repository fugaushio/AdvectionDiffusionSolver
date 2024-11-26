#include "namelist.h"

void FEMF::main_Analysis()
{
    x = 0.;
    for (int i = 0; i <= nx; i++)
    {
        u[0](i) = amplitude / 2. * exp(c / (2. * kappa) * x) * (exp(-c / (2. * kappa) * x) * erfc((x - c * t) / (2. * sqrt(kappa * t))) + exp(c / (2. * kappa) * x) * erfc((x + c * t) / (2. * sqrt(kappa * t))));
        x += dx;
    }
}