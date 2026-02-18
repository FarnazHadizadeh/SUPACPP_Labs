#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "FiniteFunctions.h"
#include "Distributions.h"

// Read 1D data points from a plain text file (one value per line).
std::vector<double> readMysteryData(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Failed to open data file: " + filename);
    }

    std::vector<double> data;
    double x = 0.0;
    while (in >> x) {
        data.push_back(x);
    }
    return data;
}

// Compute mean and standard deviation of the data.
void meanAndSigma(const std::vector<double> &data,
                  double &mean, double &sigma)
{
    if (data.empty()) {
        mean  = 0.0;
        sigma = 1.0;
        return;
    }

    double sum = 0.0;
    for (double v : data) {
        sum += v;
    }
    mean = sum / static_cast<double>(data.size());

    double sumsq = 0.0;
    for (double v : data) {
        const double d = v - mean;
        sumsq += d * d;
    }
    const double var = sumsq / static_cast<double>(data.size());
    sigma = std::sqrt(var);
    if (sigma <= 0.0) {
        sigma = 1.0;
    }
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0]
                      << " Output/Data/MysteryDataXXXXX.txt\n";
            return 1;
        }

        const std::string dataFile = argv[1];

        // -----------------------------------------------------------------
        // 1. Read the mystery data
        // -----------------------------------------------------------------
        const std::vector<double> data = readMysteryData(dataFile);
        if (data.empty()) {
            std::cerr << "No data points read from file '" << dataFile << "'.\n";
            return 1;
        }

        std::cout << "Read " << data.size()
                  << " mystery data points from " << dataFile << "\n";

        // Save a copy of the data in a simple format for plotting.
        {
            std::ofstream out("mystery_data.dat");
            for (double v : data) {
                out << v << "\n";
            }
        }

        // -----------------------------------------------------------------
        // 2. Basic statistics and sensible plotting range
        // -----------------------------------------------------------------
        const auto minmax = std::minmax_element(data.begin(), data.end());
        const double dataMin   = *minmax.first;
        const double dataMax   = *minmax.second;
        const double dataRange = dataMax - dataMin;

        double margin = 0.2 * dataRange;
        if (margin <= 0.0) {
            margin = 1.0;
        }

        const double xmin = dataMin - margin;
        const double xmax = dataMax + margin;
        const std::size_t nPoints = 2000;

        double mean  = 0.0;
        double sigma = 1.0;
        meanAndSigma(data, mean, sigma);

        std::cout << "Data range: [" << dataMin << ", " << dataMax << "]\n";
        std::cout << "Using function range: [" << xmin << ", " << xmax << "]\n";
        std::cout << "Estimated mean  = " << mean  << "\n";
        std::cout << "Estimated sigma = " << sigma << "\n\n";

        // -----------------------------------------------------------------
        // 3. Construct functions
        // -----------------------------------------------------------------
        FiniteFunctions       defaultFunc(xmin, xmax, nPoints);
        NormalFunction        normalFunc (xmin, xmax, mean, sigma, nPoints);
        CauchyLorentzFunction cauchyFunc (xmin, xmax, mean, sigma, nPoints);
        CrystalBallFunction   cbFunc     (xmin, xmax, mean, sigma,
                                          1.5, 3.0, nPoints);

        // -----------------------------------------------------------------
        // 4. Print integral / normalisation information
        // -----------------------------------------------------------------
        auto printInfo = [](const FiniteFunctions &f) {
            const double I = f.integral();
            const double N = f.normalisationFactor();
            std::cout << f.name() << "\n";
            std::cout << "  Integral over range (un-normalised) = " << I << "\n";
            std::cout << "  Normalisation factor                = " << N << "\n\n";
        };

        printInfo(defaultFunc);
        printInfo(normalFunc);
        printInfo(cauchyFunc);
        printInfo(cbFunc);

        // Write functions to file for plotting (normalised over [xmin, xmax]).
        defaultFunc.writeFunctionToFile("default_function.dat", true);
        normalFunc.writeFunctionToFile("normal_function.dat", true);
        cauchyFunc.writeFunctionToFile("cauchy_function.dat", true);
        cbFunc.writeFunctionToFile("crystalball_function.dat", true);

        std::cout << "Wrote function values to:\n"
                  << "  default_function.dat\n"
                  << "  normal_function.dat\n"
                  << "  cauchy_function.dat\n"
                  << "  crystalball_function.dat\n\n";

        // -----------------------------------------------------------------
        // 5. Metropolis sampling from the chosen distribution
        // -----------------------------------------------------------------
        std::random_device rd;
        std::mt19937 rng(rd());

        const std::size_t nSamples      = 20000;
        const std::size_t burnIn        = 1000;
        const double      proposalSigma = 0.5 * sigma; // tunable

        // Here we use the Crystal Ball as the "best guess" distribution.
        // To sample from a different function, just replace 'cbFunc' below.
        const std::vector<double> samples =
            cbFunc.metropolisSample(nSamples, proposalSigma, burnIn, rng, true);

        std::cout << "Generated " << samples.size()
                  << " Metropolis samples from " << cbFunc.name() << "\n";

        {
            std::ofstream out("metropolis_samples.dat");
            for (double v : samples) {
                out << v << "\n";
            }
        }

        std::cout << "First 10 Metropolis samples:\n";
        for (std::size_t i = 0; i < 10 && i < samples.size(); ++i) {
            std::cout << "  " << samples[i] << "\n";
        }

        std::cout << "\nWrote data files:\n"
                  << "  mystery_data.dat\n"
                  << "  metropolis_samples.dat\n"
                  << "  default_function.dat\n"
                  << "  normal_function.dat\n"
                  << "  cauchy_function.dat\n"
                  << "  crystalball_function.dat\n";

        return 0;
    } catch (const std::exception &ex) {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return 1;
    }
}
