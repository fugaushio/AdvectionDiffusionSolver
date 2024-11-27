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
    MatrixXd u_mat; // 出力用Matrix

    // vector
    VectorXd u[NUMBER];  // 現在の流速（DG法用に配列で定義．DG法以外はu[0]のみ使用）
    VectorXd u_b;        // 過去の流速
    VectorXd ue[NUMBER]; // 現在の微小要素の流速
    VectorXd ue_b;       // 過去の微小要素の流速
    // VectorXd u_a;
    // VectorXd ue_a;
    VectorXd delta_u; // ニュートン法における変分

    VectorXd gausspoint; // ガウス点={-1/sqrt(3), -1/sqrt(3)}
    double gaussweight;  // ガウスの重み=1
    double detJ;         // dx/ds=dx/2

    // matrix
    MatrixXd LHS, RHS; // 連立方程式の行列
    VectorXd R;        // 連立方程式の右辺

    // shapefunction
    double N(int num, double point);
    double dNdx(int num, double point);
    void setBC(MatrixXd &A, VectorXd &b, int flag); // 境界条件をセット．移流拡散の場合，u(0)=U0とする

    // DG
    // VectorXd ux;
    VectorXd Q[NUMBER];            // 各w_iに対応する，方程式
    MatrixXd dQdu[NUMBER][NUMBER]; // Qを各u_iで微分したもの
    double TT[NUMBER][NUMBER];     // 時間方向の形状関数をTとし，T*Tのdt区間の積分値
    double TdTdt[NUMBER][NUMBER];  // 時間方向の形状関数をTとし，T*dT/dtのdt区間の積分値

    // parameter
    double x;         // 位置
    double dx;        // 空間刻み幅
    double length;    // 全長
    int nx;           // 空間要素数
    double t;         // 時刻
    double dt;        // 時間刻み幅
    double totaltime; // 全時間
    int nt;           // 時間要素数
    int num_node_s;   // 空間要素の点の数=2
    int num_node_t;   // 時間要素の点の数=2or3

    double kappa;     // 拡散係数
    double c;         // 移流速度
    double alpha;     // Crank-nicolsonの分ける割合
    double amplitude; // 振幅=最初の左端の値
    // int BCtype;
    double tolerance; // 収束閾値
    string outputDir;
    string fileName;

    // SUPG
    double tau;

    void Hello();

    MatrixXd parameter;
    void initialize(string &filename); // パラメータの取得や，共有して使う関数の定義，行列の大きさ決めなどの行う
    void inputDomainInfo();            // 形状のパラメータを取得
    void inputBoundaryInfo();          // 境界のパラメータを取得
    void inputsolverInfo();            // 物理量のぱらめーたを取得

    // Analysis
    void main_Analysis(); // 解析解の計算
    // Explicit
    void main_Explicit(); // 陽解法（Crank-nicolson法）の計算
    // Implicit
    void main_Implicit(); // 陰解法
    // DG
    void main_DG();    // 不連続ガラーキン法
    void integral_T(); // 時間の形状関数を積分する

    // solver
    int num_threads;                                                           // 不要
    void solveLinearSystem(VectorXd &y, const MatrixXd &A, const VectorXd &b); // 連立方程式のソルバーであるが，不使用
};

#endif // namelist.h