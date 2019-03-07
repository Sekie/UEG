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
        double Lx, Ly, Lz, Volume, kMax, kF, dkx, dky, dkz, rs;
        int NumElectrons, nxMax, nyMax, nzMax;
        std::vector< std::tuple<double, int, int, int> > aOccupiedLevels;
        std::vector< std::tuple<double, int, int, int> > bOccupiedLevels;
        double EExchange, EFermi, EKinetic, n;

        void FillLevels();
        void ExciteUEG(double, double, double);
        void CalcKinetic();
        void CalcExchange();
        double CalcAnalyticalKinetic();
        double CalcAnalyticalExchange();
        double OneOverK2(int, int, bool);
        void PrintHighestOcc();
        void PrintLevel(int);
        void PrintOrbitalEnergies();

        UEG(double, double);
        UEG();
        
        double C_TF = 3.0 * M_PI * M_PI / 10.0 * pow(3.0 / M_PI, 2.0 / 3.0);
};