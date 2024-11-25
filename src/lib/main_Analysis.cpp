#include "namelist.h"

void FEMF::main_Analysis()
{
    x = 0.;
    for (int i = 0; i <= nx; i++)
    {
        if (BCtype == 0)
        {
            u[0](i) = amplitude * sin(2. * M_PI * x / length) * exp(-kappa * pow(2. * M_PI / length, 2) * t);
        }
        else if (BCtype == 1)
        {
            u[0](i) = amplitude * cos(2. * M_PI * x / length) * exp(-kappa * pow(2. * M_PI / length, 2) * t);
        }
        else
        {
            std::cout << "BCtype Error" << std::endl;
            break;
        }
        x += dx;
    }
}