#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// ======================== Base class ========================

class FiniteFunctions {
protected:
    double      m_xmin;
    double      m_xmax;
    std::size_t m_nPoints;

public:
    FiniteFunctions(double xmin, double xmax, std::size_t nPoints = 1000)
        : m_xmin(xmin), m_xmax(xmax), m_nPoints(nPoints)
    {
        if (m_xmin >= m_xmax) std::swap(m_xmin, m_xmax);
        if (m_nPoints == 0)   m_nPoints = 1000;
    }

    virtual ~FiniteFunctions() = default;

    double xmin() const { return m_xmin; }
    double xmax() const { return m_xmax; }
    std::size_t nPoints() const { return m_nPoints; }

    // Default function: f(x) = 1 / (1 + x^2)
    virtual double evaluate(double x) const {
        return 1.0 / (1.0 + x * x);
    }

    virtual std::string name() const {
        return "1 / (1 + x^2)";
    }

    // Numerical integral over [xmin, xmax] (trapezoidal rule)
    double integral(std::size_t nSamples = 0) const {
        std::size_t n = (nSamples > 0 ? nSamples : m_nPoints);
        if (n == 0) n = 1;

        double dx  = (m_xmax - m_xmin) / static_cast<double>(n);
        double sum = 0.0;

        for (std::size_t i = 0; i <= n; ++i) {
            double x = m_xmin + dx * static_cast<double>(i);
            double w = (i == 0 || i == n) ? 0.5 : 1.0;
            sum += w * evaluate(x);
        }

        return sum * dx;
    }

    // Factor to make integral = 1
    double normalisationFactor(std::size_t nSamples = 0) const {
        double I = integral(nSamples);
        if (I <= 0.0) return 1.0;
        return 1.0 / I;
    }

    // For plotting if needed: writes x f(x) to file
    void writeFunctionToFile(const std::string &filename,
                             bool normalised = true,
                             std::size_t nSamples = 0) const
    {
        std::size_t n = (nSamples > 0 ? nSamples : m_nPoints);
        if (n == 0) n = 1;

        std::ofstream out(filename.c_str());
        if (!out) return;

        double dx   = (m_xmax - m_xmin) / static_cast<double>(n);
        double norm = normalised ? normalisationFactor(nSamples) : 1.0;

        for (std::size_t i = 0; i <= n; ++i) {
            double x = m_xmin + dx * static_cast<double>(i);
            double y = evaluate(x) * norm;
            out << x << " " << y << "\n";
        }
    }

    // Metropolis sampler
    std::vector<double> metropolisSample(std::size_t nSamples,
                                         double proposalSigma,
                                         std::size_t burnIn,
                                         std::mt19937 &rng,
                                         bool normalise = true) const
    {
        std::vector<double> samples;
        samples.reserve(nSamples);
        if (nSamples == 0) return samples;

        std::uniform_real_distribution<double> uniformX(m_xmin, m_xmax);
        std::uniform_real_distribution<double> uniform01(0.0, 1.0);
        std::normal_distribution<double> proposal(0.0, proposalSigma);

        double x    = uniformX(rng);      // starting point
        double norm = normalise ? normalisationFactor() : 1.0;
        double fx   = evaluate(x) * norm;
        if (fx < 0.0) fx = 0.0;

        std::size_t totalSteps = burnIn + nSamples;

        for (std::size_t i = 0; i < totalSteps; ++i) {
            double y = x + proposal(rng);

            if (y < m_xmin || y > m_xmax) {
                // reject out of range
            } else {
                double fy = evaluate(y) * norm;
                if (fy < 0.0) fy = 0.0;

                double A;
                if (fx <= 0.0) {
                    A = 1.0;
                } else {
                    double ratio = fy / fx;
                    A = (ratio < 1.0) ? ratio : 1.0;
                }

                double T = uniform01(rng);
                if (T < A) {
                    x  = y;
                    fx = fy;
                }
            }

            if (i >= burnIn) samples.push_back(x);
        }

        return samples;
    }
};

// ======================== Distributions ========================

namespace {
    constexpr double pi = 3.14159265358979323846;
}

// --- Normal ---

class NormalFunction : public FiniteFunctions {
private:
    double m_mu;
    double m_sigma;

public:
    NormalFunction(double xmin, double xmax,
                   double mu, double sigma,
                   std::size_t nPoints = 1000)
        : FiniteFunctions(xmin, xmax, nPoints),
          m_mu(mu), m_sigma(sigma) {}

    double evaluate(double x) const override {
        const double invSqrt2Pi = 1.0 / std::sqrt(2.0 * pi);
        double z = (x - m_mu) / m_sigma;
        return invSqrt2Pi / m_sigma * std::exp(-0.5 * z * z);
    }

    std::string name() const override {
        return "Normal(mu=" + std::to_string(m_mu)
             + ", sigma=" + std::to_string(m_sigma) + ")";
    }
};

// --- Cauchy–Lorentz ---

class CauchyLorentzFunction : public FiniteFunctions {
private:
    double m_x0;
    double m_gamma;

public:
    CauchyLorentzFunction(double xmin, double xmax,
                          double x0, double gamma,
                          std::size_t nPoints = 1000)
        : FiniteFunctions(xmin, xmax, nPoints),
          m_x0(x0), m_gamma(gamma) {}

    double evaluate(double x) const override {
        double t = (x - m_x0) / m_gamma;
        return 1.0 / (pi * m_gamma * (1.0 + t * t));
    }

    std::string name() const override {
        return "Cauchy-Lorentz(x0=" + std::to_string(m_x0)
             + ", gamma=" + std::to_string(m_gamma) + ")";
    }
};

// --- Crystal Ball ---

class CrystalBallFunction : public FiniteFunctions {
private:
    double m_xbar;
    double m_sigma;
    double m_alpha;
    double m_n;

    double m_N;
    double m_A;
    double m_B;

    void updateConstants() {
        double absAlpha = std::fabs(m_alpha);
        double expTerm  = std::exp(-0.5 * absAlpha * absAlpha);

        m_A = std::pow(m_n / absAlpha, m_n) * expTerm;
        m_B = m_n / absAlpha - absAlpha;

        double C = (m_n / absAlpha) * (1.0 / (m_n - 1.0)) * expTerm;
        double D = std::sqrt(pi / 2.0) *
                   (1.0 + std::erf(absAlpha / std::sqrt(2.0)));

        m_N = 1.0 / (m_sigma * (C + D));
    }

public:
    CrystalBallFunction(double xmin, double xmax,
                        double xbar, double sigma,
                        double alpha, double n,
                        std::size_t nPoints = 1000)
        : FiniteFunctions(xmin, xmax, nPoints),
          m_xbar(xbar), m_sigma(sigma),
          m_alpha(alpha), m_n(n)
    {
        updateConstants();
    }

    double evaluate(double x) const override {
        double t = (x - m_xbar) / m_sigma;

        if (t > -m_alpha) {
            return m_N * std::exp(-0.5 * t * t);          // Gaussian core
        } else {
            return m_N * m_A * std::pow(m_B - t, -m_n);   // power tail
        }
    }

    std::string name() const override {
        return "CrystalBall(x0=" + std::to_string(m_xbar)
             + ", sigma=" + std::to_string(m_sigma)
             + ", alpha=" + std::to_string(m_alpha)
             + ", n=" + std::to_string(m_n) + ")";
    }
};

// ======================== Helpers ========================

// Generate "mystery" data: here a normal distribution
std::vector<double> generateMysteryData(int N, double mu, double sigma) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::normal_distribution<double> gauss(mu, sigma);

    std::vector<double> data;
    data.reserve(N);
    for (int i = 0; i < N; ++i) {
        data.push_back(gauss(rng));
    }
    return data;
}

// Mean and sigma of data
void meanAndSigma(const std::vector<double> &data,
                  double &mean, double &sigma)
{
    if (data.empty()) {
        mean = 0.0;
        sigma = 1.0;
        return;
    }

    double sum = 0.0;
    for (double v : data) sum += v;
    mean = sum / static_cast<double>(data.size());

    double sumsq = 0.0;
    for (double v : data) {
        double d = v - mean;
        sumsq += d * d;
    }
    double var = sumsq / static_cast<double>(data.size());
    sigma = std::sqrt(var);
}

// ======================== main ========================

int main() {
    // 1. Generate "mystery" data inside the same program
    int    N_data = 10000;
    double true_mu    = 0.5;
    double true_sigma = 1.2;

    std::vector<double> data = generateMysteryData(N_data, true_mu, true_sigma);
    std::cout << "Generated " << data.size()
              << " mystery data points from Normal("
              << true_mu << ", " << true_sigma << ")\n";

    // 2. Basic statistics and range
    auto minmax = std::minmax_element(data.begin(), data.end());
    double dataMin   = *minmax.first;
    double dataMax   = *minmax.second;
    double dataRange = dataMax - dataMin;

    double margin = 0.2 * dataRange;
    if (margin <= 0.0) margin = 1.0;

    double xmin = dataMin - margin;
    double xmax = dataMax + margin;
    std::size_t nPoints = 2000;

    double mean, sigma;
    meanAndSigma(data, mean, sigma);
    if (sigma <= 0.0) sigma = 1.0;

    std::cout << "Data range: [" << dataMin << ", " << dataMax << "]\n";
    std::cout << "Using function range: [" << xmin << ", " << xmax << "]\n";
    std::cout << "Estimated mean = " << mean
              << ", sigma = " << sigma << "\n\n";

    // 3. Construct functions
    FiniteFunctions       defaultFunc(xmin, xmax, nPoints);
    NormalFunction        normalFunc (xmin, xmax, mean, sigma, nPoints);
    CauchyLorentzFunction cauchyFunc (xmin, xmax, mean, sigma, nPoints);
    CrystalBallFunction   cbFunc     (xmin, xmax, mean, sigma,
                                      1.5, 3.0, nPoints);

    // 4. Print integral information
    auto printInfo = [&](const FiniteFunctions &f) {
        double I = f.integral();
        double N = f.normalisationFactor();
        std::cout << f.name() << "\n";
        std::cout << "  Integral over range (unnormalised) = " << I << "\n";
        std::cout << "  Normalisation factor               = " << N << "\n\n";
    };

    printInfo(defaultFunc);
    printInfo(normalFunc);
    printInfo(cauchyFunc);
    printInfo(cbFunc);

    // 5. Metropolis sampling using Crystal Ball as example
    std::random_device rd;
    std::mt19937 rng(rd());

    std::size_t nSamples     = 20000;
    std::size_t burnIn       = 1000;
    double      proposalSigma = sigma * 0.5;

    std::vector<double> samples =
        cbFunc.metropolisSample(nSamples, proposalSigma, burnIn, rng, true);

    std::cout << "Generated " << samples.size()
              << " Metropolis samples from " << cbFunc.name() << "\n";

    std::cout << "First 10 Metropolis samples:\n";
    for (std::size_t i = 0; i < 10 && i < samples.size(); ++i) {
        std::cout << samples[i] << "\n";
    }

    {
        std::ofstream outData("mystery_data.dat");
        for (double v : data) outData << v << "\n";
    }
    {
        std::ofstream outSamples("samples_cb.dat");
        for (double v : samples) outSamples << v << "\n";
    }
    cbFunc.writeFunctionToFile("crystalball_function.dat", true);

    std::cout << "\nWrote mystery_data.dat, samples_cb.dat,"
              << " and crystalball_function.dat\n";

    return 0;
}
