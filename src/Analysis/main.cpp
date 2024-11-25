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
    Femf.outputDir = "result_Analysis";
    filesystem::create_directory(Femf.outputDir);
    Femf.fileName = Femf.outputDir + "/result_" + to_string(Femf.dt);
    clock_t start = clock();
    for (int i = 0; i <= Femf.nt; i++)
    {
        Femf.main_Analysis();
        Femf.u_mat.row(i) = Femf.u[0].transpose();
        cout << "finish_T:" << Femf.t << "[s]" << endl;
        Femf.t += Femf.dt;
    }
    clock_t end = clock();
    outputCSV(Femf.u_mat, Femf.fileName);
    cout << "all finish" << endl;
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "Simulation time: " << duration << " s" << endl;
    return 0;
}
