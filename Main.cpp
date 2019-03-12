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
#include "ExcitedFunctional.h"

int main()
{
    // for (double V = 10; V < 1E6 + 1; V *= 10)
    // {
        // UEG myUEG(10.0, 1E6);
        // myUEG.CalcKinetic();
        // myUEG.CalcExchange();
        // std::cout << 1E6 << "\t" << myUEG.EKinetic / (double)myUEG.NumElectrons << "\t" << myUEG.CalcAnalyticalKinetic() / (double)myUEG.NumElectrons << "\t" << myUEG.EExchange / (double)myUEG.NumElectrons << "\t" << myUEG.CalcAnalyticalExchange() / (double)myUEG.NumElectrons << std::endl;
    // }
    double n = 10.0;
    double V = 1000.0;

    UEG myUEG(n, V);
    
    ExcitedFunctional myEF(myUEG);
    myEF.GenerateExchange();

    return 0;
}