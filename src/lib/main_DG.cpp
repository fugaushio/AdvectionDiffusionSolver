#include "namelist.h"
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

void FEMF::main_DG()
{
    u[0] = u_b;
    for (int i = 1; i < num_node_t; i++)
    {
        u[i] = Eigen::VectorXd::Zero(nx + 1);
    }

    do
    {
        // initialize
        LHS = Eigen::MatrixXd::Zero(num_node_t * (nx + 1), num_node_t * (nx + 1));
        R = Eigen::VectorXd::Zero(num_node_t * (nx + 1));
        delta_u = Eigen::VectorXd::Zero(num_node_t * (nx + 1));
        for (int i = 0; i < num_node_t; i++)
        {
            Q[i] = Eigen::VectorXd::Zero(nx + 1);
            for (int j = 0; j < num_node_t; j++)
            {
                dQdu[i][j] = Eigen::MatrixXd::Zero(nx + 1, nx + 1);
            }
        }

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < num_node_t; j++)
            {
                ue[j](0) = u[j](i);
                ue[j](1) = u[j](i + 1);
            }
            ue_b(0) = u_b(i);
            ue_b(1) = u_b(i + 1);

            double vel[NUMBER][2] = {};
            for (int i = 0; i < num_node_t; i++)
            {
                for (int j = 0; j < num_node_s; j++)
                {
                    vel[i][0] += N(j, gausspoint(0)) * ue[i](j);
                    vel[i][1] += N(j, gausspoint(1)) * ue[i](j);
                }
            }
            double vel_b[2] = {0, 0};
            for (int i = 0; i < num_node_s; i++)
            {
                vel_b[0] += N(i, gausspoint(0)) * ue_b(i);
                vel_b[1] += N(i, gausspoint(1)) * ue_b(i);
            }

            double dvdx[NUMBER][2] = {};
            for (int i = 0; i < num_node_t; i++)
            {
                for (int j = 0; j < num_node_s; j++)
                {
                    dvdx[i][0] += dNdx(j, gausspoint(0)) * ue[i](j);
                    dvdx[i][1] += dNdx(j, gausspoint(1)) * ue[i](j);
                }
            }

            for (int p = 0; p < num_node_s; p++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < num_node_t; k++)
                    {
                        for (int l = 0; l < num_node_t; l++)
                        {

                            Q[k](i + p) += TdTdt[k][l] * N(p, gausspoint(j)) * vel[l][j] * gaussweight * detJ;
                            Q[k](i + p) += c * TT[k][l] * N(p, gausspoint(j)) * dvdx[l][j] * gaussweight * detJ;
                            Q[k](i + p) += kappa * TT[k][l] * dNdx(p, gausspoint(j)) * dvdx[l][j] * gaussweight * detJ;
                        }
                    }
                    Q[0](i + p) += N(p, gausspoint(j)) * (vel[0][j] - vel_b[j]) * gaussweight * detJ;
                }
            }

            for (int p = 0; p < num_node_s; p++)
            {
                for (int q = 0; q < num_node_s; q++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        for (int k = 0; k < num_node_t; k++)
                        {
                            for (int l = 0; l < num_node_t; l++)
                            {
                                dQdu[k][l](i + p, i + q) += TdTdt[k][l] * N(p, gausspoint(j)) * N(q, gausspoint(j)) * gaussweight * detJ;
                                dQdu[k][l](i + p, i + q) += c * TT[k][l] * N(p, gausspoint(j)) * dNdx(q, gausspoint(j)) * gaussweight * detJ;
                                dQdu[k][l](i + p, i + q) += kappa * TT[k][l] * dNdx(p, gausspoint(j)) * dNdx(q, gausspoint(j)) * gaussweight * detJ;
                            }
                        }
                        dQdu[0][0](i + p, i + q) += N(p, gausspoint(j)) * N(q, gausspoint(j)) * gaussweight * detJ;
                    }
                }
            }

            // SUPG
            // tau = (pow(2. / dt, 2.) + 2. * (ue(0) * ue(0) + ue(1) * ue(1)) / dx);
            // tau = 1 / sqrt(tau);
        }
        for (int i = 0; i < num_node_t; i++)
        {
            for (int j = 0; j < num_node_t; j++)
            {
                setBC(dQdu[i][j], Q[i], 0);
                if (j != i)
                {
                    dQdu[i][j](0, 0) = 0.;
                    dQdu[i][j](0, nx) = 0.;
                }
            }
        }
        for (int k = 0; k < num_node_t; k++)
        {
            for (int i = 0; i < nx + 1; i++)
            {
                for (int l = 0; l < num_node_t; l++)
                {
                    for (int j = 0; j < nx + 1; j++)
                    {
                        LHS(i + k * (nx + 1), j + l * (nx + 1)) = dQdu[k][l](i, j);
                    }
                }
                R(i + k * (nx + 1)) = -Q[k](i);
            }
        }
        delta_u = LHS.colPivHouseholderQr().solve(R);
        // solveLinearSystem(delta_u, LHS, R);
        for (int k = 0; k < num_node_t; k++)
        {
            for (int i = 0; i < nx + 1; i++)
            {
                u[k](i) += delta_u(i + k * (nx + 1));
            }
        }
        // std::cout << "Norm:" << delta_u.norm() << endl;
        if (delta_u.norm() > 1e1)
        {
            exit(1);
            std::cout << "Error" << std::endl;
        }
    } while (delta_u.norm() > tolerance);
    for (int i = 1; i < num_node_t; i++)
    {
        u[0] += u[i];
    }
}
