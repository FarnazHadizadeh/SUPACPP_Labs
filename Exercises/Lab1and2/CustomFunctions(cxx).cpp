#include "CustomFunctions.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>

// ---------------- Reading ----------------
/*
XYData ReadXYFile(const std::string& filename) {
    XYData data;
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Error: could not open data file '" << filename << "'.\n";
        return data;
    }

    double x, y;
    while (in >> x >> y) {   // standard two-column read
        data.emplace_back(x, y);
    }
    return data;
} */

XYData ReadXYFile(const std::string& filename) {
	XYData data;
	std::ifstream in(filename);
	if (!in) {
		std::cerr << "Error: could not open data file '" << filename << "'.\n";
		return data;
	}

	double x, y;
	while (true) {
		// Read x; if this fails, we are done.
		if (!(in >> x)) break;

		// If the next non-extracted char is a comma, consume it.
		int next = in.peek();
		if (next == ',') {
			in.get();  // remove the comma
		}

		// Now read y; if this fails, we are done (incomplete pair).
		if (!(in >> y)) break;

		data.emplace_back(x, y);
	}

	return data;
}

/*
std::vector<double> ReadErrorFile(const std::string& filename) {
    std::vector<double> errs;
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Warning: could not open error file '" << filename
                  << "'. Chi2 will assume sigma=1 for all points.\n";
        return errs; // empty signals fallback
    }

    // Robust: allow 1 or 2 columns; take last value in each line.
    std::string line;
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        double val, last = 0.0;
        bool any = false;
        while (ss >> val) { last = val; any = true; }
        if (any) errs.push_back(last);
    }
    return errs;
} */

std::vector<double> ReadErrorFile(const std::string& filename) {
    std::vector<double> errs;
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Warning: could not open error file '" << filename
                  << "'. Chi2 will assume sigma=1 for all points.\n";
        return errs; // empty signals fallback
    }

    std::string line;
    while (std::getline(in, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');

        std::stringstream ss(line);
        double val, last = 0.0;
        bool any = false;
        while (ss >> val) {
            last = val;
            any = true;
        }
        if (any) errs.push_back(last);
    }
    return errs;
}


// ---------------- Printing ----------------

void PrintFirstN(const XYData& data, int N) {
	if (data.empty()) {
		std::cerr << "Warning: no data loaded.\n";
		return;
	}
	if (N <= 0) {
		std::cerr << "Warning: N must be positive.\n";
		return;
	}

	int toPrint = N;
	if (static_cast<size_t>(N) > data.size()) {
		std::cerr << "Warning: N larger than data size (" << data.size()
		          << "). Printing first 5 lines only.\n";
		toPrint = std::min<int>(5, static_cast<int>(data.size()));
	}

	for (int i = 0; i < toPrint; ++i) {
		std::cout << i << ": x=" << data[i].first
		          << ", y=" << data[i].second << "\n";
	}
}

void PrintVector(const std::vector<double>& v) {
	for (size_t i = 0; i < v.size(); ++i) {
		std::cout << i << ": " << v[i] << "\n";
	}
}

void PrintVector(const std::vector<long double>& v) {
	for (size_t i = 0; i < v.size(); ++i) {
		std::cout << i << ": " << v[i] << "\n";
	}
}

// ---------------- Magnitudes ----------------

std::vector<double> ComputeMagnitudes(const XYData& data) {
	std::vector<double> mags;
	mags.reserve(data.size());

	for (const auto& pt : data) {
		double x = pt.first, y = pt.second;
		mags.push_back(std::sqrt(x*x + y*y));
	}
	return mags;
}

// ---------------- Least squares + chi2 ----------------

std::pair<double,double> FitLineLeastSquares(const XYData& data) {
	const size_t N = data.size();
	if (N < 2) {
		std::cerr << "Error: need at least 2 points for a line fit.\n";
		return {0.0, 0.0};
	}

	long double sumx = 0.0L, sumy = 0.0L, sumxx = 0.0L, sumxy = 0.0L;
	for (const auto& pt : data) {
		long double x = pt.first, y = pt.second;
		sumx  += x;
		sumy  += y;
		sumxx += x*x;
		sumxy += x*y;
	}

	long double denom = static_cast<long double>(N)*sumxx - sumx*sumx;
	if (std::fabs(denom) < 1e-18L) {
		std::cerr << "Error: denominator nearly zero; cannot fit line.\n";
		return {0.0, 0.0};
	}

	long double p = (static_cast<long double>(N)*sumxy - sumx*sumy) / denom;
	long double q = (sumxx*sumy - sumx*sumxy) / denom;

	return {static_cast<double>(p), static_cast<double>(q)};
}

double Chi2NDOF(const XYData& data,
                const std::vector<double>& errors,
                double p, double q) {

	const size_t N = data.size();
	if (N == 0) return 0.0;

	long double chi2 = 0.0L;

	for (size_t i = 0; i < N; ++i) {
		double x = data[i].first;
		double y_obs = data[i].second;
		double y_exp = p*x + q;

		double sigma = 1.0;
		if (!errors.empty() && i < errors.size()) sigma = errors[i];
		if (sigma == 0.0) {
			std::cerr << "Warning: sigma=0 at i=" << i << ", skipping.\n";
			continue;
		}

		long double resid = y_obs - y_exp;
		chi2 += (resid*resid) / (sigma*sigma);
	}

	long double ndof = static_cast<long double>(N) - 2.0L;
	if (ndof <= 0) return static_cast<double>(chi2);

	return static_cast<double>(chi2 / ndof);
}

// ---------------- Recursive power + x^y ----------------

long double PowerRecursive(long double base, long long exp) {
	if (exp == 0) return 1.0L;
	if (exp < 0)  return 1.0L / PowerRecursive(base, -exp);

	// exponentiation by squaring (recursive, no loops)
	if (exp % 2 == 0) {
		long double half = PowerRecursive(base, exp/2);
		return half * half;
	} else {
		return base * PowerRecursive(base, exp-1);
	}
}

std::vector<long double> ComputeXPowY(const XYData& data) {
	std::vector<long double> out;
	out.reserve(data.size());

	std::transform(data.begin(), data.end(), std::back_inserter(out),
	[](const auto& pt) {
		long double x = pt.first;
		long long y = static_cast<long long>(std::llround(pt.second));
		return PowerRecursive(x, y);
	});

	return out;
}

// ---------------- Saving results ----------------

void SaveToFile(const std::string& filename, const XYData& data) {
	std::ofstream out(filename);
	if (!out) {
		std::cerr << "Error: cannot write to '" << filename << "'.\n";
		return;
	}
	for (const auto& pt : data) {
		out << pt.first << " " << pt.second << "\n";
	}
}

void SaveToFile(const std::string& filename, const std::vector<double>& data) {
	std::ofstream out(filename);
	if (!out) {
		std::cerr << "Error: cannot write to '" << filename << "'.\n";
		return;
	}
	for (double v : data) out << v << "\n";
}

void SaveToFile(const std::string& filename, const std::vector<long double>& data) {
	std::ofstream out(filename);
	if (!out) {
		std::cerr << "Error: cannot write to '" << filename << "'.\n";
		return;
	}
	for (long double v : data) out << v << "\n";
}

void SaveToFile(const std::string& filename, const std::string& text) {
	std::ofstream out(filename);
	if (!out) {
		std::cerr << "Error: cannot write to '" << filename << "'.\n";
		return;
	}
	out << text << "\n";
}
