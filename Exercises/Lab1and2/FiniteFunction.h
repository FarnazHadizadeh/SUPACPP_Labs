#ifndef FINITEFUNCTIONS_H
#define FINITEFUNCTIONS_H

#include <cstddef>
#include <random>
#include <string>
#include <vector>

class FiniteFunctions {
protected:
    double      m_xmin;
    double      m_xmax;
    std::size_t m_nPoints;

public:
    /// Construct a finite function defined on [xmin, xmax] with a default
    /// number of points used for numerical integration / plotting.
    explicit FiniteFunctions(double xmin, double xmax, std::size_t nPoints = 1000);

    virtual ~FiniteFunctions() = default;

    double xmin() const;
    double xmax() const;
    std::size_t nPoints() const;

    /// Default function: f(x) = 1 / (1 + x^2)
    virtual double evaluate(double x) const;

    /// Human–readable name of the function.
    virtual std::string name() const;

    /// Numerical integral over [xmin, xmax] using the trapezoidal rule.
    double integral(std::size_t nSamples = 0) const;

    /// Multiplicative factor such that ∫ f(x) * normalisationFactor dx = 1.
    double normalisationFactor(std::size_t nSamples = 0) const;

    /// Write (x, f(x)) pairs to file for plotting.
    /// If normalised == true the function values are multiplied by normalisationFactor().
    void writeFunctionToFile(const std::string &filename,
                             bool normalised = true,
                             std::size_t nSamples = 0) const;

    /// Sample from the (optionally normalised) function using the Metropolis algorithm.
    /// The proposal is a normal distribution with standard deviation 'proposalSigma'.
    /// 'burnIn' steps are discarded. The returned vector has size nSamples.
    std::vector<double> metropolisSample(std::size_t nSamples,
                                         double proposalSigma,
                                         std::size_t burnIn,
                                         std::mt19937 &rng,
                                         bool normalise = true) const;
};

#endif // FINITEFUNCTIONS_H
