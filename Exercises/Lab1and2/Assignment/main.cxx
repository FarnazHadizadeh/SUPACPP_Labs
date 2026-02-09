#include <iostream>
#include "FiniteFunctions.h"
#include "CustomDistributions.h"

int main() {
    int npts = 10000;

    FiniteFunctions base(-10, 10, npts);
    NormalDistribution normal(-10, 10, npts, 0.0, 1.0);
    CauchyDistribution cauchy(-10, 10, npts, 0.0, 1.0);

    std::cout << "Base Integral: " << base.Integrate() << std::endl;
    std::cout << "Normal Integral: " << normal.Integrate() << std::endl;
    std::cout << "Cauchy Integral: " << cauchy.Integrate() << std::endl;

    return 0;
}
