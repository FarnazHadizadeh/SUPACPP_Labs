#include "FiniteFunctions.h"
#include <cmath>

FiniteFunctions::FiniteFunctions(double xmin, double xmax, int npts)
    : xMin(xmin), xMax(xmax), nPoints(npts) {}

FiniteFunctions::~FiniteFunctions() {}

double FiniteFunctions::Evaluate(double x) const {
    return 1.0 / (1.0 + x * x);
}

double FiniteFunctions::Integrate() const {
    double dx = (xMax - xMin) / nPoints;
    double sum = 0.0;

    for (int i = 0; i < nPoints; i++) {
        double x = xMin + i * dx;
        sum += Evaluate(x);
    }
    return sum * dx;
}

std::vector<double> FiniteFunctions::GetXValues() const {
    std::vector<double> xvals;
    double dx = (xMax - xMin) / nPoints;
    for (int i = 0; i < nPoints; i++)
        xvals.push_back(xMin + i * dx);
    return xvals;
}

std::vector<double> FiniteFunctions::GetYValues() const {
    std::vector<double> yvals;
    for (double x : GetXValues())
        yvals.push_back(Evaluate(x));
    return yvals;
}

double FiniteFunctions::GetXMin() const { return xMin; }
double FiniteFunctions::GetXMax() const { return xMax; }
