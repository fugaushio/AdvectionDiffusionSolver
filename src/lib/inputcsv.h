#pragma once
#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

Eigen::MatrixXi readCSV_i(const std::string &filename)
{
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file");
    }

    std::vector<std::vector<int>> data;
    std::string line;

    // Skip the header line
    if (!std::getline(file, line))
    {
        throw std::runtime_error("No data found in file");
    }

    while (std::getline(file, line))
    {
        std::vector<int> row;
        std::stringstream lineStream(line);
        std::string cell;

        while (std::getline(lineStream, cell, ','))
        {
            row.push_back(std::stod(cell));
        }
        data.push_back(row);
    }

    file.close();

    if (data.empty())
    {
        throw std::runtime_error("No data found in file");
    }

    // Determine the number of rows and columns in the CSV
    size_t rows = data.size();
    size_t cols = data[0].size();

    Eigen::MatrixXi mat(rows, cols);

    // Copy data to Eigen matrix
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            mat(i, j) = data[i][j];
        }
    }

    return mat;
}

Eigen::MatrixXd readCSV_d(const std::string &filename)
{
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file");
    }

    std::vector<std::vector<double>> data;
    std::string line;

    // Skip the header line
    if (!std::getline(file, line))
    {
        throw std::runtime_error("No data found in file");
    }

    while (std::getline(file, line))
    {
        std::vector<double> row;
        std::stringstream lineStream(line);
        std::string cell;

        while (std::getline(lineStream, cell, ','))
        {
            row.push_back(std::stod(cell));
        }
        data.push_back(row);
    }

    file.close();

    if (data.empty())
    {
        throw std::runtime_error("No data found in file");
    }

    // Determine the number of rows and columns in the CSV
    size_t rows = data.size();
    size_t cols = data[0].size();

    Eigen::MatrixXd mat(rows, cols);

    // Copy data to Eigen matrix
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            mat(i, j) = data[i][j];
        }
    }

    return mat;
}