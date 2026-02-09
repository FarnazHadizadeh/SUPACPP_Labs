#ifndef FINITEFUNCTIONS_H
#define FINITEFUNCTIONS_H

#include <vector>
#include <string>

class FiniteFunctions {
protected:
    double xMin, xMax;
    int nPoints;

public:
    FiniteFunctions(double xmin, double xmax, int npts);
    virtual ~FiniteFunctions();

    virtual double Evaluate(double x) const;
    double Integrate() const;

    std::vector<double> GetXValues() const;
    std::vector<double> GetYValues() const;

    double GetXMin() const;
    double GetXMax() const;
};

#endif
