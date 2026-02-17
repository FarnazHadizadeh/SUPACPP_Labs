#include "CustomFunctions.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>

// Compute magnitudes sqrt(x^2 + y^2)
std::vector<double> computeMagnitudes(const std::vector<double>& x_values,
                                      const std::vector<double>& y_values)
{
    std::vector<double> magnitudes;
    magnitudes.reserve(x_values.size());

    for (std::size_t i = 0; i < x_values.size(); ++i) {
        double mag = std::sqrt(x_values[i] * x_values[i] +
                               y_values[i] * y_values[i]);
        magnitudes.push_back(mag);
    }

    return magnitudes;
}

// Linear regression y = m x + c
void linearFit(const std::vector<double>& x_values,
               const std::vector<double>& y_values,
               double& m, double& c)
{
    double sumX  = 0.0;
    double sumY  = 0.0;
    double sumXY = 0.0;
    double sumX2 = 0.0;

    const int n = static_cast<int>(x_values.size());

    for (int i = 0; i < n; ++i) {
        sumX  += x_values[i];
        sumY  += y_values[i];
        sumXY += x_values[i] * y_values[i];
        sumX2 += x_values[i] * x_values[i];
    }

    m = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    c = (sumY - m * sumX) / n;
}

// Print contents of a vector
void print(const std::vector<double>& values)
{
    for (double v : values) {
        std::cout << v << std::endl;
    }
}

// Read error values from file
std::vector<double> readErrors(const std::string& filename)
{
    std::vector<double> errors;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: could not open error file " << filename << "\n";
        return errors;
    }

    double e;
    while (file >> e) {
        errors.push_back(e);
    }

    return errors;
}

// Compute chi-squared per degree of freedom
double computeChi2NDOF(const std::vector<double>& x_values,
                       const std::vector<double>& y_values,
                       const std::vector<double>& errors,
                       double m, double c)
{
    std::size_t n = std::min({ x_values.size(), y_values.size(), errors.size() });

    if (n <= 2) {
        std::cerr << "Not enough points to compute chi2/NDOF.\n";
        return 0.0;
    }

    double chi2 = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        double y_model = m * x_values[i] + c;
        double diff    = y_values[i] - y_model;
        chi2 += (diff * diff) / (errors[i] * errors[i]);
    }

    int ndof = static_cast<int>(n) - 2; // two fit parameters: m and c
    if (ndof <= 0) {
        std::cerr << "Non-positive NDOF.\n";
        return 0.0;
    }

    return chi2 / ndof;
}
