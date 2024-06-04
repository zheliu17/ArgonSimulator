#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include "header.h"

int main(int argc, char *argv[])
{
    std::vector<std::function<void()>> functionList = {
        test_neighbor_list,
        test_neighbor_force,
        test_pbc,
        benchmark_neighbor_freq_nopbc,
        benchmark_neighbor_freq_pbcon,
        test_omp,
        NVT_results};

    if (argc > 1)
    {
        size_t functionIndex = std::stoi(argv[1]);

        if (functionIndex >= 0 && functionIndex < functionList.size())
        {
            functionList[functionIndex]();
        }
        else
        {
            std::cerr << "Invalid function index\n";
        }
    }
    else
    {
        std::cerr << "Please provide a function index as an argument\n";
    }


    return 0;
}
