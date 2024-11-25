#include "namelist.h"
#include <vector>
#include "Eigen/Dense"

using namespace Eigen;


void FEMF::main_Implicit()
{
    do
    {
        LHS = Eigen::MatrixXd::Zero(nx + 1, nx + 1);
        R = Eigen::VectorXd::Zero(nx + 1);
        delta_u = Eigen::VectorXd::Zero(nx + 1);
        for (int i = 0; i < nx; i++)
        {
            ue[0](0) = u[0](i);
            ue[0](1) = u[0](i + 1);
            ue_b(0) = u_b(i);
            ue_b(1) = u_b(i + 1);

            double vel[2] = {0, 0};
            for (int i = 0; i < num_node_s; i++)
            {
                vel[0] += N(i, gausspoint(0)) * ue[0](i);
                vel[1] += N(i, gausspoint(1)) * ue[0](i);
            }
            double vel_b[2] = {0, 0};
            for (int i = 0; i < num_node_s; i++)
            {
                vel_b[0] += N(i, gausspoint(0)) * ue_b(i);
                vel_b[1] += N(i, gausspoint(1)) * ue_b(i);
            }
            double dvdx[2] = {0, 0};
            for (int i = 0; i < num_node_s; i++)
            {
                dvdx[0] += dNdx(i, gausspoint(0)) * ue[0](i);
                dvdx[1] += dNdx(i, gausspoint(1)) * ue[0](i);
            }

            for (int p = 0; p < num_node_s; p++)
            {
                for (int j = 0; j < 2; j++)
                {
                    R(i + p) -= N(p, gausspoint(j)) * (vel[j] - vel_b[j]) * gaussweight * detJ / dt;
                    R(i + p) -= kappa * dNdx(p, gausspoint(j)) * dvdx[j] * gaussweight * detJ;
                }
            }

            for (int p = 0; p < num_node_s; p++)
            {
                for (int q = 0; q < num_node_s; q++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        LHS(i + p, i + q) += N(p, gausspoint(j)) * N(q, gausspoint(j)) * gaussweight * detJ / dt;
                        LHS(i + p, i + q) += kappa * dNdx(p, gausspoint(j)) * dNdx(q, gausspoint(j)) * gaussweight * detJ;
                    }
                }
            }

            // SUPG
            // tau = (pow(2. / dt, 2.) + 2. * (ue(0) * ue(0) + ue(1) * ue(1)) / dx);
            // tau = 1 / sqrt(tau);
        }

        setBC(LHS, R, BCtype);
        delta_u = LHS.colPivHouseholderQr().solve(R);
        // solveLinearSystem(delta_u, LHS, R);
        u[0] += delta_u;
        // std::cout << "Norm:" << delta_u.norm() << endl;
    } while (delta_u.norm() > tolerance);
}
