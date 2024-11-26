#include "namelist.h"

void FEMF::main_Analysis()
{
    x = 0.;
    for (int i = 0; i <= nx; i++)
    {
        if (BCtype == 0)
        {
            double A = -(t + 1.) / (4. * kappa * t);
            double B = (length * length / 4. * (-t * t + t + 1) - length * t * (x - c * t) + pow(x - c * t, 2) / t) / pow(t + 1., 2);
            std::cout<<"A:"<<A<<" B:"<<B<<std::endl;
            u[0](i) = amplitude / sqrt(4. * M_PI * kappa * t) * exp(A * B) * sqrt(M_PI / A);
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