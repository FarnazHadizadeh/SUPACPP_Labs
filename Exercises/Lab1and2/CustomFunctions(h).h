#ifndef CUSTOMFUNCTIONS_H
#define CUSTOMFUNCTIONS_H

#include <string>
#include <vector>
#include <utility>

// Type alias for clarity
using XYData = std::vector<std::pair<double,double>>;

// Reading
XYData ReadXYFile(const std::string& filename);
std::vector<double> ReadErrorFile(const std::string& filename);

// Printing
void PrintFirstN(const XYData& data, int N);
void PrintVector(const std::vector<double>& v);
void PrintVector(const std::vector<long double>& v);

// Magnitudes
std::vector<double> ComputeMagnitudes(const XYData& data);

// Least squares + chi2
std::pair<double,double> FitLineLeastSquares(const XYData& data);
double Chi2NDOF(const XYData& data,
                const std::vector<double>& errors,
                double p, double q);

// Recursive power + x^y
long double PowerRecursive(long double base, long long exp);
std::vector<long double> ComputeXPowY(const XYData& data);

// Saving (overloads)
void SaveToFile(const std::string& filename, const XYData& data);
void SaveToFile(const std::string& filename, const std::vector<double>& data);
void SaveToFile(const std::string& filename, const std::vector<long double>& data);
void SaveToFile(const std::string& filename, const std::string& text);

#endif
