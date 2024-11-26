#ifndef _NAMELIST_H_
#define _NAMELIST_H_
#define NUMBER 3

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <omp.h>
#include <array>

// F#include "TextParser.h"
//  #include "petsc_cpp.h"
#include <Eigen/Core>
#include <Eigen/Dense>
// #include <Eigen/LU>
// #include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

class FEMF
{
public:
    // output u Matrix
    MatrixXd u_mat;

    // vector
    VectorXd u[NUMBER];
    VectorXd u_b;
    VectorXd ue[NUMBER];
    VectorXd ue_b;
    // VectorXd u_a;
    // VectorXd ue_a;
    VectorXd delta_u;

    VectorXd gausspoint;
    double gaussweight;
    double detJ;

    // matrix
    MatrixXd LHS, RHS;
    VectorXd R;

    // shapefunction
    double N(int num, double point);
    double dNdx(int num, double point);
    void setBC(MatrixXd &A, VectorXd &b, int flag);

    // DG
    // VectorXd ux;
    VectorXd Q[NUMBER];
    MatrixXd dQdu[NUMBER][NUMBER];
    double TT[NUMBER][NUMBER];
    double TdTdt[NUMBER][NUMBER];

    // parameter
    double x;
    double dx;
    double length;
    int nx;
    double t;
    double dt;
    double totaltime;
    int nt;
    int num_node_s;
    int num_node_t;

    double kappa;
    double c;
    double alpha;
    double amplitude;
    // int BCtype;
    double tolerance;
    string outputDir;
    string fileName;

    // SUPG
    double tau;

    void Hello();

    MatrixXd parameter;
    void initialize(string &filename);
    void inputDomainInfo();
    void inputBoundaryInfo();
    void inputsolverInfo();

    // Analysis
    void main_Analysis();
    // Explicit
    void main_Explicit();
    // Implicit
    void main_Implicit();
    // DG
    void main_DG();
    void integral_T();

    // solver
    int num_threads;
    void solveLinearSystem(VectorXd &y, const MatrixXd &A, const VectorXd &b);
};

#endif // namelist.h