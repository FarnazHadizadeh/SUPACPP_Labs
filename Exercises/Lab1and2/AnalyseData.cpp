#include "CustomFunctions.h"

#include <iostream>
#include <string>

int main() {
    std::cout << "=== AnalyseData ===\n";

    // Allow user to override filename, but provide a safe default.
    std::string dataFile = "input2D_float.txt";
    std::cout << "Enter data filename (default: input2D_float.txt): ";
    std::string tmp;
    std::getline(std::cin >> std::ws, tmp);
    if (!tmp.empty()) dataFile = tmp;

    XYData data = ReadXYFile(dataFile);
    if (data.empty()) {
        std::cerr << "No data loaded. Exiting.\n";
        return 1;
    }
    std::cout << "Loaded " << data.size() << " points.\n\n";


    bool running = true;
    while (running) {
        std::cout
            << "\nChoose an option:\n"
            << "  1) Print N lines of (x,y)\n"
            << "  2) Compute magnitudes |r|\n"
            << "  3) Least-squares line fit + chi2/NDOF\n"
            << "  4) Compute x^y (recursive, y rounded)\n"
            << "  5) Exit\n"
            << "Selection: ";

        int choice;
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::cin.ignore(100000, '\n');
            std::cerr << "Invalid input.\n";
            continue;
        }

        switch (choice) {
            case 1: {
                int N;
                std::cout << "Enter N: ";
                std::cin >> N;
                PrintFirstN(data, N);
                SaveToFile("output_printN.txt", data);
                std::cout << "Saved full data to output_printN.txt\n";
                break;
            }
            case 2: {
                auto mags = ComputeMagnitudes(data);
                std::cout << "Magnitudes:\n";
                PrintVector(mags);
                SaveToFile("output_magnitudes.txt", mags);
                std::cout << "Saved to output_magnitudes.txt\n";
                break;
            }
            case 3: {
                auto [p, q] = FitLineLeastSquares(data);
                auto errs = ReadErrorFile("error2D_float.txt");
                double chi2ndof = Chi2NDOF(data, errs, p, q);

                std::string fitText =
                    "Fit: y = p x + q\n"
                    "p = " + std::to_string(p) + "\n" +
                    "q = " + std::to_string(q) + "\n" +
                    "chi2/NDOF = " + std::to_string(chi2ndof) + "\n";

                std::cout << fitText;
                SaveToFile("output_fit.txt", fitText);
                std::cout << "Saved to output_fit.txt\n";
                break;
            }
            case 4: {
                auto xpowy = ComputeXPowY(data);
                std::cout << "x^y results:\n";
                PrintVector(xpowy);
                SaveToFile("output_xpowy.txt", xpowy);
                std::cout << "Saved to output_xpowy.txt\n";
                break;
            }
            case 5:
                running = false;
                break;

            default:
                std::cerr << "Unknown option.\n";
        }

        if (running) {
            std::cout << "\nRun another task? (y/n): ";
            char again;
            std::cin >> again;
            if (again != 'y' && again != 'Y') running = false;
        }
    }

    std::cout << "Exiting AnalyseData.\n";
    return 0;
}
