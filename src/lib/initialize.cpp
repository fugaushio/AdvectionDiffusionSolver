#include "namelist.h"
#include "inputcsv.h"
#include <filesystem>

using namespace std;
using namespace Eigen;

void FEMF::Hello()
{
    std::cout << "HelloWorld" << std::endl;
    return;
}

void FEMF::initialize(string &filename)
{
    parameter = readCSV_d(filename);
    inputDomainInfo();
    inputBoundaryInfo();
    inputsolverInfo();

    nx = length / dx;
    nt = totaltime / dt;
    for (int i = 0; i < num_node_t; i++)
    {
        u[i] = Eigen::VectorXd::Zero(nx + 1);
        ue[i] = Eigen::VectorXd::Zero(num_node_s);
    }

    u_b = Eigen::VectorXd::Zero(nx + 1);
    // u_a = Eigen::VectorXd::Zero(nx + 1);

    ue_b = Eigen::VectorXd::Zero(num_node_s);
    // ue_a = Eigen::VectorXd::Zero(num_node_s);
    x = 0.;
    t = 0.;
    u_mat = Eigen::MatrixXd::Zero(nt + 1, nx + 1);

    // shapefunction0

    integral_T();

    gausspoint = Eigen::VectorXd::Zero(2);
    gausspoint(0) = -1 / sqrt(3.);
    gausspoint(1) = 1 / sqrt(3);
    gaussweight = 1.;
    detJ = dx / 2.;

    u_b(0) = amplitude;
    for (int i = 1; i <= nx; i++)
    {
        u_b(i) = 0.;
        x += dx;
    }
}

void FEMF::inputDomainInfo()
{
    dx = parameter(0, 0);
    length = parameter(0, 1);
    dt = parameter(0, 2);
    totaltime = parameter(0, 3);
}
void FEMF::inputBoundaryInfo()
{
    // BCtype = static_cast<int>(parameter(0, 4));
    amplitude = parameter(0, 4);
}
void FEMF::inputsolverInfo()
{
    kappa = parameter(0, 5);
    c = parameter(0, 6);
    tolerance = parameter(0, 7);
    num_node_s = parameter(0, 8);
    num_node_t = parameter(0, 9);
    alpha = parameter(0, 10);
}

void FEMF::integral_T()
{
    if (num_node_t == 2)
    {
        TT[0][0] = 1. * dt;
        TT[0][1] = 1. / 2. * dt;
        TT[1][0] = 1. / 2. * dt;
        TT[1][1] = 1. / 3. * dt;
        TdTdt[0][0] = 0.;
        TdTdt[0][1] = 1.;
        TdTdt[1][0] = 0.;
        TdTdt[1][1] = 1. / 2.;
    }
    else if (num_node_t == 3)
    {
        TT[0][0] = 1. * dt;
        TT[0][1] = 5. / 6. * dt;
        TT[0][2] = 1. / 6. * dt;
        TT[1][0] = 5. / 6. * dt;
        TT[1][1] = 4. / 5. * dt;
        TT[1][2] = 1. / 5. * dt;
        TT[2][0] = 1. / 6. * dt;
        TT[2][1] = 1. / 5. * dt;
        TT[2][2] = 2. / 15 * dt;
        TdTdt[0][0] = 0.;
        TdTdt[0][1] = 1.;
        TdTdt[0][2] = 1.;
        TdTdt[1][0] = 0.;
        TdTdt[1][1] = 1. / 2.;
        TdTdt[1][2] = 7. / 6.;
        TdTdt[2][0] = 0.;
        TdTdt[2][1] = -1. / 6.;
        TdTdt[2][2] = 1. / 2.;
    }
}

void FEMF::solveLinearSystem(VectorXd &y, const MatrixXd &A, const VectorXd &b)
{
    // SparseMatrix<double> sparseA = A.sparseView();
    // SparseLU<SparseMatrix<double>> solver;
    // solver.compute(sparseA);
    // if (solver.info() == Success)
    //    y = solver.solve(b);
    // else
    y = A.fullPivLu().solve(b);
}

double FEMF::N(int num, double point)
{
    double N0, N1;
    N0 = (1. - point) / 2.;
    N1 = (1. + point) / 2.;
    if (num == 0)
        return N0;
    else if (num == 1)
        return N1;
    else
        return 0;
}

double FEMF::dNdx(int num, double point)
{

    double dNds0, dNds1;
    dNds0 = -1. / 2.;
    dNds1 = 1. / 2.;
    if (num == 0)
        return dNds0 / detJ;
    else if (num == 1)
        return dNds1 / detJ;
    else
        return 0;
}

void FEMF::setBC(MatrixXd &A, VectorXd &b, int flag)
{
    int n = A.rows() - 1;
    if (flag == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            A(0, i) = 0.;
            // A(n, i) = 0.;
        }
        A(0, 0) = 1.;
        // A(0, n) = -1.;

        b(0) = 0.;
        // b(n) = 0.;
    }

    if (flag == 1)
    {
        for (int i = 0; i < n + 1; i++)
        {
            A(0, i) = 0.;
            // A(n, i) = 0.;
        }
        A(0, 0) = 1.;
        // A(0, n) = -1.;

        b(0) = amplitude;
        // b(n) = 0.;
    }
}