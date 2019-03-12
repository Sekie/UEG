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
    std::ofstream Output("xvsk.txt");
    for (int n = 10; n < 100; n++)
    {
        for (int i = 0; i < 10; i++)
        {
            myUEG.RandomExciteUEG(n, n, n);
            myUEG.CalcKinetic();
            myUEG.CalcExchange();
            Output << myUEG.EKinetic / myUEG.Volume << "\t" << myUEG.EExchange / (double)myUEG.NumElectrons << std::endl;
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
