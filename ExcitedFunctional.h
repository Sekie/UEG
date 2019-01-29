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

class ExcitedFunctional
{
    public:
        double rs;
        

        double CalculateExchange();
    
    private:
        double C_TF = 0.3 * M_PI * M_PI * pow(3 / M_PI, 2.0 / 3.0);
};