#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <tuple>
#include <fstream>
#include <map>
#include <stdlib.h> 
#include <algorithm> // std::sort
#include <iomanip>
#include <queue>
#include <unsupported/Eigen/CXX11/Tensor>

#include "UEG.h"

void UEG::FillLevels()
{
    int NumPairs = NumElectrons / 2;
    std::vector< std::tuple<double, int, int, int> > UnsortedLevels;
    int SearchN = 100;
    for (int nx = 0; nx < SearchN; nx++)
    {
        for (int ny = 0; ny < SearchN; ny++)
        {
            for (int nz = 0; nz < SearchN; nz++)
            {
                double E = 2.0 * M_PI * M_PI * ((nx / Lx) * (nx / Lx) + (ny / Ly) * (ny / Ly) + (nz / Lz) * (nz / Lz));
                std::tuple<double, int, int, int> tmpTuple = std::make_tuple(E, nx, ny, nz);
                UnsortedLevels.push_back(tmpTuple);
            }
        }
    }
    std::sort(UnsortedLevels.begin(), UnsortedLevels.end());
    for (int i = 0; i < NumPairs; i++)
    {
        aOccupiedLevels.push_back(UnsortedLevels[i]);
        bOccupiedLevels.push_back(UnsortedLevels[i]);
    }
    if (NumElectrons % 2 != 0)
    {
        aOccupiedLevels.push_back(UnsortedLevels[NumPairs]);
    }
}

void UEG::CalcKinetic()
{
    EKinetic = 0.0;
    for (int i = 0; i < aOccupiedLevels.size(); i++)
    {
        EKinetic += std::get<0>(aOccupiedLevels[i]);
    }
    for (int i = 0; i < bOccupiedLevels.size(); i++)
    {
        EKinetic += std::get<0>(bOccupiedLevels[i]);
    }
}

double UEG::OneOverK2(int i, int j, bool SpinAlpha)
{
    double dnx, dny, dnz;
    if (SpinAlpha)
    {
        dnx = std::get<1>(aOccupiedLevels[i]) - std::get<1>(aOccupiedLevels[j]);
        dny = std::get<2>(aOccupiedLevels[i]) - std::get<2>(aOccupiedLevels[j]);
        dnz = std::get<3>(aOccupiedLevels[i]) - std::get<3>(aOccupiedLevels[j]);
    }
    else
    {
        dnx = std::get<1>(bOccupiedLevels[i]) - std::get<1>(bOccupiedLevels[j]);
        dny = std::get<2>(bOccupiedLevels[i]) - std::get<2>(bOccupiedLevels[j]);
        dnz = std::get<3>(bOccupiedLevels[i]) - std::get<3>(bOccupiedLevels[j]);
    }

    double KInv = 4.0 * M_PI * M_PI * (dnx * dnx / (Lx * Lx) + dny * dny / (Ly * Ly) + dnz * dnz / (Lz * Lz));
    KInv = 1.0 / KInv;

    
    return KInv;
}

void UEG::CalcExchange()
{
    EExchange = 0.0;
    int NumPairs = NumElectrons / 2;
    for (int i = 0; i < NumPairs; i++)
    {
        for (int j = i + 1; j < NumPairs; j++)
        {
            EExchange += OneOverK2(i, j, true) + OneOverK2(i, j, false);
        }
    }

    if (NumElectrons % 2 != 0)
    {
        for (int j = 0; j < NumPairs; j++)
        {
            EExchange += OneOverK2(NumPairs, j, true);
        }
    }

    EExchange *= -4.0 * M_PI / Volume;
}

UEG::UEG(int NumElec, double V)
{
    NumElectrons = NumElec;
    Volume = V;
    Lx = Ly = Lz = cbrt(Volume);
    n = (double)NumElectrons / Volume;
}
