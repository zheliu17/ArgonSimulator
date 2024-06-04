#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

using std::vector;

int nearest_cubic(int number)
{
    int power = static_cast<int>(std::ceil(std::log(number) / std::log(3)));
    int candidate = power * power * power;
    if (candidate <= number)
    {
        while (power * power * power < number)
        {
            power += 1;
        }
    }
    else
    {
        while (power * power * power > number)
        {
            power -= 1;
        }
        power += 1;
    }
    return power;
}


void sortVectorsByFirst(
    vector<int> &firstVector,
    vector<vector<double>> &secondVector,
    vector<int> &thirdVector)
{
    // Create an index vector to maintain the original indices
    std::vector<std::size_t> indices(firstVector.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Define a custom comparator that compares elements in the first vector
    auto comparator = [&firstVector](std::size_t a, std::size_t b)
    {
        return firstVector[a] < firstVector[b];
    };

    // Sort indices based on the values in the first vector
    std::sort(indices.begin(), indices.end(), comparator);

    // Rearrange both vectors based on the sorted indices
    vector<int> sortedFirstVector(firstVector.size());
    vector<vector<double>> sortedSecondVector(secondVector.size());
    vector<int> sortedThirdVector(thirdVector.size());

    for (std::size_t i = 0; i < indices.size(); ++i)
    {
        sortedFirstVector[i] = firstVector[indices[i]];
        sortedSecondVector[i] = secondVector[indices[i]];
        sortedThirdVector[i] = thirdVector[indices[i]];
    }
    firstVector = sortedFirstVector;
    secondVector = sortedSecondVector;
    thirdVector = sortedThirdVector;
}