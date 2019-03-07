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

#include "ExcitedFunctional.h"
#include "UEG.h"

double ExcitedFunctional::CalculateExchange()
{
    double rho = 4 * M_PI / (3 * rs * rs * rs);

}

void ExcitedFunctional::GenerateExchange()
{
    std::ofstream Output("out.txt");
    for (int i = 0; i < 100; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            double kg = (double)i / 100.0 * myUEG.kF;
            double kx = ((double)j / 100.0 + 1.0) * myUEG.kF;
            myUEG.ExciteUEG(kg, 0.01, kx);
            myUEG.CalcKinetic();
            myUEG.CalcExchange();
            Output << myUEG.EKinetic << "\t" << myUEG.EExchange << std::endl;
        }
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