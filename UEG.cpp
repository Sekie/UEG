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
    double k2Max = kF * kF;
    for (int nx = -nxMax; nx < nxMax + 1; nx++)
    {
        for (int ny = -nyMax; ny < nyMax + 1; ny++)
        {
            for (int nz = -nzMax; nz < nzMax + 1; nz++)
            {
                double k2 = nx * nx * dkx * dkx + ny * ny * dky * dky + nz * nz * dkz * dkz;
                if (k2 < k2Max)
                {
                    double E = 0.5 * k2;
                    std::tuple<double, int, int, int> tmpTuple = std::make_tuple(E, nx, ny, nz);
                    aOccupiedLevels.push_back(tmpTuple);
                    bOccupiedLevels.push_back(tmpTuple);
                }
            }
        }
    }
    std::sort(aOccupiedLevels.begin(), aOccupiedLevels.end());
    std::sort(bOccupiedLevels.begin(), bOccupiedLevels.end());
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
    if (fabs(KInv) < 1E-3)
    {
        std::cout << "KInv is 0 for " << i << " " << j << " " << SpinAlpha << std::endl;
    }
    KInv = 1.0 / KInv;
    
    return KInv;
}

void UEG::CalcExchange()
{
    EExchange = 0.0;
    for (int i = 0; i < aOccupiedLevels.size(); i++)
    {
        for (int j = i + 1; j < aOccupiedLevels.size(); j++)
        {
            EExchange += OneOverK2(i, j, true);
        }
    }
    // for (int i = 0; i < bOccupiedLevels.size(); i++)
    // {
    //     for (int j = i + 1; j < bOccupiedLevels.size(); j++)
    //     {
    //         EExchange += OneOverK2(i, j, false);
    //     }
    // }

    EExchange *= 2.0 * -4.0 * M_PI / Volume;
}

double UEG::CalcAnalyticalKinetic()
{
    double KE = 0.6 * std::get<0>(aOccupiedLevels[aOccupiedLevels.size() - 1]) * NumElectrons;
    return KE;
}

double UEG::CalcAnalyticalExchange()
{
    double X = -Volume * pow(kF, 4) / (4.0 * M_PI * M_PI * M_PI);
    return X;
}

void UEG::PrintHighestOcc()
{
    std::cout << "The highest occupied level is:" << std::endl;
    std::cout << std::get<1>(aOccupiedLevels[aOccupiedLevels.size() - 1]) << "\t" << std::get<2>(aOccupiedLevels[aOccupiedLevels.size() - 1]) << "\t" << std::get<3>(aOccupiedLevels[aOccupiedLevels.size() - 1]) << std::endl;
}

void UEG::PrintLevel(int Level)
{
    std::cout << "The " << Level << " level is:" << std::endl;
    std::cout << std::get<1>(aOccupiedLevels[Level]) << "\t" << std::get<2>(aOccupiedLevels[Level]) << "\t" << std::get<3>(aOccupiedLevels[Level]) << std::endl;
}

void UEG::PrintOrbitalEnergies()
{
    std::cout << "The orbital energies are:" << std::endl;
    for (int i = 0; i < NumElectrons / 2; i++)
    {
        std::cout << std::get<0>(aOccupiedLevels[i]) << std::endl;
    }
}

void UEG::ExciteUEG(double kg, double dk, double kx)
{
    int NumOcc = aOccupiedLevels.size();
    aOccupiedLevels.clear();
    bOccupiedLevels.clear();
    std::vector< std::tuple<double, int, int, int> > ExcitedLevels;
    double k2Max = kF * kF;
    for (int nx = -4 * nxMax; nx < 4 * nxMax + 1; nx++)
    {
        for (int ny = -4 * nyMax; ny < 4 * nyMax + 1; ny++)
        {
            for (int nz = -4 * nzMax; nz < 4 * nzMax + 1; nz++)
            {
                double k2 = nx * nx * dkx * dkx + ny * ny * dky * dky + nz * nz * dkz * dkz;
                if ((k2 < k2Max && k2 < kg * kg && k2 > (kg + dk) * (kg + dk)))
                {
                    double E = 0.5 * k2;
                    std::tuple<double, int, int, int> tmpTuple = std::make_tuple(E, nx, ny, nz);
                    aOccupiedLevels.push_back(tmpTuple);
                    bOccupiedLevels.push_back(tmpTuple);
                }

                if (k2 > kx * kx)
                {
                    double E = 0.5 * k2;
                    std::tuple<double, int, int, int> tmpTuple = std::make_tuple(E, nx, ny, nz);
                    ExcitedLevels.push_back(tmpTuple);
                }
            }
        }
    }
    std::sort(aOccupiedLevels.begin(), aOccupiedLevels.end());
    std::sort(bOccupiedLevels.begin(), bOccupiedLevels.end());
    std::sort(ExcitedLevels.begin(), ExcitedLevels.end());
    
    int MissingOcc = NumOcc - aOccupiedLevels.size();
    for (int i = 0; i < MissingOcc; i++)
    {
        aOccupiedLevels.push_back(ExcitedLevels[i]);
        bOccupiedLevels.push_back(ExcitedLevels[i]);
    }
}

void UEG::RandomExciteUEG(int MaxNx, int MaxNy, int MaxNz)
{
    int NumOcc = aOccupiedLevels.size();
    aOccupiedLevels.clear();
    bOccupiedLevels.clear();

    for (int i = 0; i < NumOcc; i++)
    {
        int nx = (rand() % (2 * MaxNx)) - MaxNx;
        int ny = (rand() % (2 * MaxNy)) - MaxNy;
        int nz = (rand() % (2 * MaxNz)) - MaxNz;
        bool Repeat = false;
        for (int j = 0; j < aOccupiedLevels.size(); j++)
        {
            if (nx == std::get<1>(aOccupiedLevels[j]) && ny == std::get<2>(aOccupiedLevels[j]) && nz == std::get<3>(aOccupiedLevels[j]))
            {
                Repeat = true;
                break;
            }
        }
        if (Repeat)
        {
            i--;
            continue;
        }
        // std::cout << nx << "\t" << ny << "\t" << nz << std::endl;
        double k2 = nx * nx * dkx * dkx + ny * ny * dky * dky + nz * nz * dkz * dkz;
        double E = k2 / 2.0;
        std::tuple<double, int, int, int> tmpTuple = std::make_tuple(E, nx, ny, nz);
        aOccupiedLevels.push_back(tmpTuple);
        bOccupiedLevels.push_back(tmpTuple);
    }
}

UEG::UEG(double p, double V)
{
    NumElectrons = p * V;
    Volume = V;
    Lx = Ly = Lz = cbrt(Volume);
    n = p;
    rs = pow(3.0 / (4.0 * M_PI * p), 1.0 / 3.0);

    kF = pow(3.0 * M_PI * M_PI * p, 1.0 / 3.0);
    dkx = dky = dkz = 2 * M_PI / Lx;
    nxMax = nyMax = nzMax = ceil(kF / dkx);

    FillLevels();
}

UEG::UEG()
{
    // Just a default constructor
}
