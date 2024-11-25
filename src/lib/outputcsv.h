#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

void outputCSV(const Eigen::MatrixXd &matrix, const std::string &filename)
{
    std::ofstream file(filename + ".csv"); // ファイル名に ".csv" を追加

    if (file.is_open())
    {
        for (int i = 0; i < matrix.rows(); ++i)
        {
            for (int j = 0; j < matrix.cols(); ++j)
            {
                file << matrix(i, j);
                if (j < matrix.cols() - 1) // 列の区切りにカンマを追加
                    file << ",";
            }
            file << "\n"; // 行の区切りに改行を追加
        }
        file.close();
    }
    else
    {
        std::cerr << "Could not open file\n";
    }
}