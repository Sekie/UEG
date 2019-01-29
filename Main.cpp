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
    UEG myUEG(1E2, 1E4);
    myUEG.FillLevels();
    myUEG.CalcKinetic();
    std::cout << 0.6 * std::get<0>(myUEG.aOccupiedLevels[myUEG.aOccupiedLevels.size() - 1]) << std::endl;
    std::cout << std::get<1>(myUEG.aOccupiedLevels[myUEG.aOccupiedLevels.size() - 1]) << "\t" << std::get<2>(myUEG.aOccupiedLevels[myUEG.aOccupiedLevels.size() - 1]) << "\t" << std::get<3>(myUEG.aOccupiedLevels[myUEG.aOccupiedLevels.size() - 1]) << std::endl;
    std::cout << myUEG.EKinetic / (double)myUEG.NumElectrons << std::endl;

    myUEG.CalcExchange();
    std::cout << myUEG.EExchange << std::endl;
}