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

int main()
{
    double n = 10.0;
    double V = 10000.0;

    double rs = pow(3.0 / (4.0 * M_PI * n), 1.0 / 3.0);
    std::cout << "rs = " << rs << std::endl;
    UEG myUEG(n, V);
    myUEG.PrintHighestOcc();
    myUEG.CalcKinetic();
    std::cout << 0.6 * std::get<0>(myUEG.aOccupiedLevels[myUEG.aOccupiedLevels.size() - 1]) << std::endl;
    std::cout << myUEG.EKinetic / (double)myUEG.NumElectrons << std::endl;

    myUEG.CalcExchange();
    std::cout << -V * pow(std::get<0>(myUEG.aOccupiedLevels[myUEG.aOccupiedLevels.size() - 1]), 4) / (4.0 * M_PI * M_PI * M_PI * (double)myUEG.NumElectrons) << std::endl;
    // std::cout << - pow(std::get<0>(myUEG.aOccupiedLevels[myUEG.aOccupiedLevels.size() - 1]), 4) << std::endl;
    std::cout << myUEG.EExchange << std::endl;

    std::cout << "E_HF / N = " << 2.21 / (rs * rs) << " - " << 0.916 / rs << " = " << 2.21 / (rs * rs) - 0.916 / rs << std::endl;
    std::cout << "E_EX / N = " << 2.21 / (rs * rs) - 0.916 / rs - 0.6 * std::get<0>(myUEG.aOccupiedLevels[myUEG.aOccupiedLevels.size() - 1]) << std::endl;
}