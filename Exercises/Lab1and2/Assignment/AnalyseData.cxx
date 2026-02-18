#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "CustomFunctions.h"

int main()
{
    std::vector<double> x_values, y_values;
    std::vector<double> errors;

    int choice = 0;

    while (choice != 5) {
        std::cout << "\n===== MENU =====\n";
        std::cout << "1. Load data file\n";
        std::cout << "2. Print first N lines\n";
        std::cout << "3. Compute magnitudes\n";
        std::cout << "4. Linear fit (least squares)\n";
        std::cout << "5. Exit\n";
        std::cout << "Enter choice: ";
        std::cin >> choice;

        // 1) Load data file
        if (choice == 1) {
            std::string filename;
            std::cout << "Enter filename (e.g. Data/MysteryData20000.txt): ";
            std::cin >> filename;

            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cout << "Error: could not open file " << filename << "\n";
                continue;
            }

            x_values.clear();
            y_values.clear();

            double x, y;
            while (file >> x >> y) {
                x_values.push_back(x);
                y_values.push_back(y);
            }

            std::cout << "File loaded: " << x_values.size() << " rows.\n";

            // Load error values for chi2/NDOF
            errors = readErrors("Data/error2D_float.txt");
            if (errors.empty()) {
                std::cout << "Warning: error file could not be loaded or is empty.\n";
            } else {
                std::cout << "Error file loaded: " << errors.size()
                          << " values.\n";
            }
        }

        // 2) Print first N lines
        else if (choice == 2) {
            if (x_values.empty() || y_values.empty()) {
                std::cout << "No data loaded. Please choose option 1 first.\n";
                continue;
            }

            int N;
            std::cout << "Enter N: ";
            std::cin >> N;

            for (int i = 0; i < N && i < static_cast<int>(x_values.size()); ++i) {
                std::cout << x_values[i] << "   " << y_values[i] << std::endl;
            }
        }

        // 3) Compute magnitudes
        else if (choice == 3) {
            if (x_values.empty() || y_values.empty()) {
                std::cout << "No data loaded. Please choose option 1 first.\n";
                continue;
            }

            std::vector<double> magnitudes =
                computeMagnitudes(x_values, y_values);

            std::cout << "Magnitudes:\n";
            print(magnitudes);
        }

        // 4) Linear fit and chi2/NDOF
        else if (choice == 4) {
            if (x_values.empty() || y_values.empty()) {
                std::cout << "No data loaded. Please choose option 1 first.\n";
                continue;
            }

            double m, c;
            linearFit(x_values, y_values, m, c);
            std::cout << "Linear fit result: y = " << m << " * x + " << c << "\n";

            if (!errors.empty()) {
                double chi2_ndof =
                    computeChi2NDOF(x_values, y_values, errors, m, c);
                std::cout << "Chi2 / NDOF = " << chi2_ndof << "\n";
            } else {
                std::cout
                    << "Chi2 / NDOF could not be computed (no error data loaded).\n";
            }
        }

        // 5) Exit
        else if (choice == 5) {
            std::cout << "Exiting.\n";
        }

        // Invalid choice
        else {
            std::cout << "Invalid choice, please try again.\n";
        }
    }

    return 0;
}
