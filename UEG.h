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
        std::vector< std::tuple<double, int, int, int> > aVirtualLevels;
        std::vector< std::tuple<double, int, int, int> > bVirtualLevels;
        std::vector< std::vector< double > > InvK2;
        double EExchange, EFermi, EKinetic, VarExchange, n;

        std::vector<int> OccIdx;
        std::vector<int> VirIdx;

        void FillLevels();
        void ExciteUEG(double, double, double);
        void RandomExciteUEG(int, int, int);
        void GetVirtual();
        void StoreInvK2();
        void CalcKinetic();
        void CalcExchange();
        void CalcExchangeVar();
        void CalcExchangeVarFromStorage();
        double CalcAnalyticalKinetic();
        double CalcAnalyticalExchange();
        double OneOverK2(int, int, bool);
        double OneOverK2OccVir(int, int, bool);
        void PrintHighestOcc();
        void PrintLevel(int);
        void PrintOrbitalEnergies();

        void SetNMax(int);

        UEG(double, double);
        UEG();
        
        double C_TF = 3.0 * M_PI * M_PI / 10.0 * pow(3.0 / M_PI, 2.0 / 3.0);
};