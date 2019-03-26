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
#include <numeric>

#include "ExcitedFunctional.h"
#include "UEG.h"

double CalcVariance(std::vector<double> Vec)
{
    double Sum = std::accumulate(Vec.begin(), Vec.end(), 0.0);
    double Mean = Sum / Vec.size();
    double Variance = 0.0;
    for (int i = 0; i < Vec.size(); i++)
    {
        Variance += ((Vec[i] - Mean) * (Vec[i] - Mean));
    }
    Variance /= (Vec.size() - 1);
    return Variance;
}

double ExcitedFunctional::CalculateExchange()
{
    double rho = 4 * M_PI / (3 * rs * rs * rs);

}

void ExcitedFunctional::GenerateExchange()
{
    std::ofstream Output("xvsk.txt");
    std::ofstream OutputVar("xvarvsk.txt");
    for (int n = 10; n < 30; n++)
    {
        myUEG.SetNMax(n);

        myUEG.RandomExciteUEG(n, n, n);
        myUEG.GetVirtual();
        myUEG.CalcKinetic();
        myUEG.CalcExchange();
        Output << myUEG.EKinetic / myUEG.Volume << "\t" << myUEG.EExchange / (double)myUEG.NumElectrons << std::endl;

        myUEG.CalcExchangeVar();
        OutputVar << myUEG.EKinetic / myUEG.Volume << "\t" << myUEG.VarExchange << std::endl;
    }
}

void ExcitedFunctional::GenerateExchangeByP()
{
    double V = myUEG.Volume;
    double kg = 0.5 * myUEG.kF;
    double kx = 1.2 * myUEG.kF;
    std::ofstream pOutput("exchange_by_p.txt");
    for (int i = 1; i < 10001; i++)
    {
        double n = (double)i / 100.0 * 10.0;
        UEG UEGit(n, V);
        UEGit.ExciteUEG(kg, 0.01, kx);
        UEGit.CalcExchange();
        pOutput << n << "\t" << UEGit.EExchange << std::endl;
    }
}

ExcitedFunctional::ExcitedFunctional(UEG UEGObj)
{
    myUEG = UEGObj;
}
