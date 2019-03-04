#pragma once

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

class UEG
{
    public:
        double Lx, Ly, Lz, Volume, kMax, kF, dkx, dky, dkz;
        int NumElectrons, nxMax, nyMax, nzMax;
        std::vector< std::tuple<double, int, int, int> > aOccupiedLevels;
        std::vector< std::tuple<double, int, int, int> > bOccupiedLevels;
        double EExchange, EFermi, EKinetic, n;

        void FillLevels();
        void CalcKinetic();
        void CalcExchange();
        double OneOverK2(int, int, bool);
        void PrintHighestOcc();
        void PrintLevel(int);
        void PrintOrbitalEnergies();

        UEG(double, double);
};