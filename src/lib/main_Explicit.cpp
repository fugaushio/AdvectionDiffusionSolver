#include "namelist.h"
#include <vector>
#include "Eigen/Dense"

using namespace Eigen;

void FEMF::main_Explicit()
{
    LHS = Eigen::MatrixXd::Zero(nx + 1, nx + 1);
    RHS = Eigen::MatrixXd::Zero(nx + 1, nx + 1);
    R = Eigen::VectorXd::Zero(nx + 1);

    for (int i = 0; i < nx; i++)
    {
        for (int p = 0; p < num_node_s; p++)
        {
            for (int q = 0; q < num_node_s; q++)
            {
                for (int j = 0; j < 2; j++)
                {
                    LHS(i + p, i + q) += N(p, gausspoint(j)) * N(q, gausspoint(j)) * gaussweight * detJ / dt;
                    LHS(i + p, i + q) += kappa * alpha * dNdx(p, gausspoint(j)) * dNdx(q, gausspoint(j)) * gaussweight * detJ;
                    RHS(i + p, i + q) += N(p, gausspoint(j)) * N(q, gausspoint(j)) * gaussweight * detJ / dt;
                    LHS(i + p, i + q) -= c * N(p, gausspoint(j)) * dNdx(q, gausspoint(j)) * gaussweight * detJ;
                    RHS(i + p, i + q) -= kappa * (1. - alpha) * dNdx(p, gausspoint(j)) * dNdx(q, gausspoint(j)) * gaussweight * detJ;
                }
            }
        }
    }

    R = RHS * u_b;
    setBC(LHS, R, 0);
    u[0] = LHS.colPivHouseholderQr().solve(R);
    // solveLinearSystem(u, LHS, R);
}
