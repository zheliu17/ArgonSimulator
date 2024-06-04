#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <algorithm>
#include "header.h"

using std::vector;

void clean_neighbor_list(vector<Molecule> &molList)
{
    for (Molecule &mol : molList)
    {
        mol.neighbor_list = {};
        mol.is_pbc = {};
        mol.image_position = {};
    }
}

void full_neighbor_list(vector<Molecule> &molList)
/*Full neighbor list, one mol can see all rest mols*/
{
    int num_mol = static_cast<int>(molList.size());
    vector<int> baseVector(num_mol);
    iota(baseVector.begin(), baseVector.end(), 0);
    for (int i = 0; i < num_mol; i++)
    {
        vector<int> currentVector(baseVector);
        currentVector.erase(currentVector.begin() + i);
        molList[i].neighbor_list = currentVector;
    }
}

void half_neighbor_list(vector<Molecule> &molList)
/*Half neighbor list to avoid double counted force and energy*/
/*doesn't work with the pbc condition*/
{
    size_t num_mol = molList.size();
    for (size_t i = 0; i < num_mol; i += 1)
    {
        vector<int> currentVector(num_mol - i - 1);
        iota(currentVector.begin(), currentVector.end(), i + 1);
        molList[i].neighbor_list = currentVector;
    }
}

GridSearch::GridSearch(vector<Molecule> &molList,
                       double &R,
                       bool is_pbc_on,
                       double box_size)
    : molList(molList),
      extend_Rc(R),
      is_pbc_on(is_pbc_on),
      box_size(box_size),
      Rg(R * 0.5),
      numBins(static_cast<int>(std::ceil(box_size / Rg))),
      binIndices(
          numBins, vector<vector<vector<int>>>(
                       numBins, vector<vector<int>>(numBins))),
      r2(extend_Rc * extend_Rc), inv_Rg(1 / Rg) {}

void GridSearch::kernel()
// This function will first clean the neighbor list, then build a new binIndices
{
    clean_neighbor_list(molList);
    clean_binIndices();

    for (size_t i = 0; i < molList.size(); ++i)
    {
        int binIndexX = getBinIndex(molList[i].atoms[0].position[0], inv_Rg, numBins);
        int binIndexY = getBinIndex(molList[i].atoms[0].position[1], inv_Rg, numBins);
        int binIndexZ = getBinIndex(molList[i].atoms[0].position[2], inv_Rg, numBins);
        binIndices[binIndexX][binIndexY][binIndexZ].push_back(i);
    }

    vector<double> temp_vec = {};

#pragma omp parallel for collapse(3)
    for (int i = 0; i < numBins; ++i)
    {
        for (int j = 0; j < numBins; ++j)
        {
            for (int k = 0; k < numBins; ++k)
            {
                vector<int> current_box = binIndices[i][j][k];
                // below are 62 boxes of 125 total box
                // we only search half of 125 to avoid the double counting
                processNeighborsForI(i, j, k, current_box);
                processNeighborsForJ(i, j, k, current_box);
                processNeighborsForK(i, j, k, current_box);
                for (size_t ii = 0; ii < current_box.size(); ++ii)
                {
                    // all particles in the center box should belong to the neighbor list
                    // to avoid the double counting, we only do half
                    for (size_t jj = 0; jj < ii; ++jj)
                    {
                        molList[current_box[ii]].neighbor_list.push_back(current_box[jj]);
                        if (is_pbc_on)
                        {
                            molList[current_box[ii]].image_position.push_back(temp_vec);
                            molList[current_box[ii]].is_pbc.push_back(1);
                        }
                    }
                    if (is_pbc_on)
                    {
                        // sort the neighbor list for better cache read
                        sortVectorsByFirst(
                            molList[current_box[ii]].neighbor_list,
                            molList[current_box[ii]].image_position,
                            molList[current_box[ii]].is_pbc);
                    }
                    else
                    {
                        // sort the neighbor list for better cache read
                        std::sort(molList[current_box[ii]].neighbor_list.begin(),
                                  molList[current_box[ii]].neighbor_list.end());
                    }
                }
            }
        }
    }
}

void GridSearch::clean_binIndices()
{
    for (auto &i : binIndices)
    {
        for (auto &j : i)
        {
            for (auto &k : j)
            {
                k.clear();
            }
        }
    }
}

int GridSearch::getBinIndex(double &position, double inv_Rn, int numBins)
// Build binIndex, and reset positions if pbc
// We first reset position then do neighbor search
{
    int binIndex = static_cast<int>(position * inv_Rn);
    if (binIndex < 0)
    {
        binIndex = 0;
    }
    else if (binIndex >= numBins)
    {
        binIndex = numBins - 1;
    }
    if (is_pbc_on)
    {
        if (position < 0)
        {
            position += box_size;
        }
        else if (position > box_size)
        {
            position -= box_size;
        }
    }
    return binIndex;
}

void GridSearch::processNeighbor(int &neighbor_i, int &neighbor_j, int &neighbor_k, const vector<int> &current_box)
// Function to find the correct neighbor_box for neighbor box index (i, j, k)
{
    double off_set_x = 0;
    double off_set_y = 0;
    double off_set_z = 0;
    vector<int> neighbor_box = {};
    if (is_pbc_on)
    {
        // if pbc, update out-of-bound-ed box index to the correct one
        // and update off_set
        int temp_i = getPBCIndex(neighbor_i, off_set_x);
        int temp_j = getPBCIndex(neighbor_j, off_set_y);
        int temp_k = getPBCIndex(neighbor_k, off_set_z);
        neighbor_box = binIndices[temp_i][temp_j][temp_k];
    }
    // if not pbc, only do in-bound box
    else if (neighbor_i >= 0 && neighbor_i < numBins &&
             neighbor_j >= 0 && neighbor_j < numBins &&
             neighbor_k >= 0 && neighbor_k < numBins)
    {
        neighbor_box = binIndices[neighbor_i][neighbor_j][neighbor_k];
    }
    // if not pbc and out-of-bound this neighbor_box is empty
    search_neighbor_box(current_box, neighbor_box, off_set_x, off_set_y, off_set_z);
}

int GridSearch::getPBCIndex(const int &binIndex, double &off_set)
{
    int new_idx = binIndex;
    if (binIndex < 0)
    {
        new_idx = binIndex + numBins;
        off_set = box_size;
    }
    else if (binIndex >= numBins)
    {
        new_idx = binIndex - numBins;
        off_set = -box_size;
    }
    return new_idx;
}

// The following 3 functions search total 62 of 125 boxes

void GridSearch::processNeighborsForI(int i, int j, int k, const vector<int> &current_box)
// 50 of 125 neighbor boxes
{
    for (int neighbor_i = i + 1; neighbor_i < i + 3; ++neighbor_i)
    {
        for (int neighbor_j = j - 2; neighbor_j < j + 3; ++neighbor_j)
        {
            for (int neighbor_k = k - 2; neighbor_k < k + 3; ++neighbor_k)
            {
                processNeighbor(neighbor_i, neighbor_j, neighbor_k, current_box);
            }
        }
    }
}

void GridSearch::processNeighborsForJ(int i, int j, int k, const vector<int> &current_box)
// 10 of 125 neighbor boxes
{
    int neighbor_i = i;
    for (int neighbor_j = j + 1; neighbor_j < j + 3; ++neighbor_j)
    {
        for (int neighbor_k = k - 2; neighbor_k < k + 3; ++neighbor_k)
        {
            processNeighbor(neighbor_i, neighbor_j, neighbor_k, current_box);
        }
    }
}

void GridSearch::processNeighborsForK(int i, int j, int k, const vector<int> &current_box)
// 2 of 125 neighbor boxes
{
    int neighbor_i = i;
    int neighbor_j = j;
    for (int neighbor_k = k + 1; neighbor_k < k + 3; ++neighbor_k)
    {
        processNeighbor(neighbor_i, neighbor_j, neighbor_k, current_box);
    }
}

void GridSearch::search_neighbor_box(const vector<int> &current_box, const vector<int> &neighbor_box,
                                     const double &off_set_x, const double &off_set_y, const double &off_set_z)
{
    for (int i : current_box)
    {
        const vector<double> &pos0 = molList[i].atoms[0].position;
        double x0 = pos0[0] + off_set_x;
        double y0 = pos0[1] + off_set_y;
        double z0 = pos0[2] + off_set_z;
        for (int j : neighbor_box)
        {
            const vector<double> &pos1 = molList[j].atoms[0].position;
            double x1 = pos1[0];
            double y1 = pos1[1];
            double z1 = pos1[2];
            if ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0) < r2)
            {
                molList[i].neighbor_list.push_back(j);
                if (is_pbc_on)
                {
                    vector<double> temp_vec;
                    if (off_set_x != 0 || off_set_y != 0 || off_set_z != 0)
                    // this neighbor molecule is a image
                    {
                        temp_vec = {-off_set_x, -off_set_y, -off_set_z};
                        molList[i].is_pbc.push_back(2);
                    }
                    // this neighbor molecule is regular
                    else
                    {
                        // push empty vector
                        // ensure image_position and is_pbc have the same size
                        temp_vec = {};
                        molList[i].is_pbc.push_back(1);
                    }
                    molList[i].image_position.push_back(temp_vec);
                }
            }
        }
    }
}