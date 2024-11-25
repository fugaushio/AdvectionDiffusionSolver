#include <iostream>
#include "namelist.h"
#include "outputcsv.h"
#include <filesystem>
#include <time.h>
using namespace std;

int main(int argc, const char **argv)
{
    string inputfile = argv[1];
    FEMF Femf;
    Femf.initialize(inputfile);
    Femf.outputDir = "result_Explicit";
    filesystem::create_directory(Femf.outputDir);
    Femf.fileName = Femf.outputDir + "/result_" + to_string(Femf.dt);

    clock_t start = clock();
    Femf.u_mat.row(0) = Femf.u_b.transpose();
    cout << "finish_T:" << Femf.t << "[s]" << endl;
    Femf.t += Femf.dt;
    for (int i = 1; i <= Femf.nt; i++)
    {
        Femf.main_Explicit();
        Femf.u_mat.row(i) = Femf.u[0].transpose();
        cout << "finish_T:" << Femf.t << "[s]" << endl;
        Femf.t += Femf.dt;
        Femf.u_b = Femf.u[0];
    }
    clock_t end = clock();
    outputCSV(Femf.u_mat, Femf.fileName);
    cout << "all finish" << endl;
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "Simulation time: " << duration << " s" << endl;
    return 0;
}
