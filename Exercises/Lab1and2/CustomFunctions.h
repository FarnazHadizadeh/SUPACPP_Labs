#ifndef CUSTOMFUNCTIONS_H
#define CUSTOMFUNCTIONS_H

#include <vector>
#include <string>

// Compute magnitudes sqrt(x^2 + y^2)
std::vector<double> computeMagnitudes(const std::vector<double>& x_values,
                                      const std::vector<double>& y_values);

// Linear regression y = m x + c
void linearFit(const std::vector<double>& x_values,
               const std::vector<double>& y_values,
               double& m, double& c);

// Print contents of a vector, one value per line
void print(const std::vector<double>& values);

// Read error values from file
std::vector<double> readErrors(const std::string& filename);

// Compute chi-squared per degree of freedom
double computeChi2NDOF(const std::vector<double>& x_values,
                       const std::vector<double>& y_values,
                       const std::vector<double>& errors,
                       double m, double c);

#endif

