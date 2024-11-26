#include "namelist.h"

void FEMF::main_Analysis()
{
    x = 0.;
    for (int i = 0; i <= nx; i++)
    {

        u[0](i) = amplitude / 2. * sin(2. * M_PI * (x - c * t) - 4. * M_PI * M_PI * kappa * t);
        x += dx;
    }
}