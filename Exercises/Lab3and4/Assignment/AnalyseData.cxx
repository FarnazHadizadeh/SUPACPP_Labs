
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

// ==================== TASK FUNCTIONS ====================

// Task 1: Print first N lines
void task1(const vector<double>& X, const vector<double>& Y) {
    ofstream outfile("output_part1.txt");
    int N;
    cout << "Enter N: ";
    cin >> N;

    outfile << "First " << N << " lines of data:\n";
    for (int i = 0; i < N && i < X.size(); i++) {
        outfile << X[i] << " " << Y[i] << "\n";
    }

    outfile.close();
    cout << "Results written to output_part1.txt\n";
}

// Task 2: Calculate magnitudes
void task2(const vector<double>& X, const vector<double>& Y) {
    ofstream outfile("output_part2.txt");

    outfile << "Magnitudes:\n";
    for (int i = 0; i < X.size(); i++) {
        double mag = sqrt(X[i]*X[i] + Y[i]*Y[i]);
        outfile << "(" << X[i] << ", " << Y[i] << ") -> " << mag << "\n";
    }

    outfile.close();
    cout << "Results written to output_part2.txt\n";
}

// Task 3: Least squares line fit
void task3(const vector<double>& X, const vector<double>& Y) {
    ofstream outfile("output_part3.txt");

    int n = X.size();
    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;

    for (int i = 0; i < n; i++) {
        sumX += X[i];
        sumY += Y[i];
        sumXY += X[i] * Y[i];
        sumXX += X[i] * X[i];
    }

    double p = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
    double q = (sumXX * sumY - sumXY * sumX) / (n * sumXX - sumX * sumX);

    outfile << "Least squares line:\n";
    outfile << "y = " << p << "x + " << q << "\n";

    double chi2 = 0;
    for (int i = 0; i < n; i++) {
        double Ei = p * X[i] + q;
        chi2 += (Y[i] - Ei) * (Y[i] - Ei);
    }

    outfile << "Chi^2 = " << chi2 << "\n";
    outfile.close();
    cout << "Results written to output_part3.txt\n";
}

// Task 4: Calculate x^y without loops/pow
void task4(const vector<double>& X, const vector<double>& Y) {
    ofstream outfile("output_part4.txt");

    outfile << "x^y results:\n";
    for (int i = 0; i < X.size(); i++) {
        double result = exp(Y[i] * log(X[i])); // using exp(y * ln(x))
        outfile << X[i] << "^" << Y[i] << " = " << result << "\n";
    }

    outfile.close();
    cout << "Results written to output_part4.txt\n";
}

// Task 5: Template function for magnitude
template <typename T>
T magnitude(T x, T y) {
    return sqrt(x*x + y*y);
}

void task5() {
    ofstream outfile("output_part5.txt");

    double x, y;
    cout << "Enter x: ";
    cin >> x;
    cout << "Enter y: ";
    cin >> y;

    double mag = magnitude(x, y);
    outfile << "User-entered vector magnitude:\n";
    outfile << "Magnitude = " << mag << "\n";

    outfile.close();
    cout << "Results written to output_part5.txt\n";
}

// ==================== MAIN ====================
int main() {
    string filename = "data.txt";
    ifstream infile(filename);

    if (!infile) {
        cout << "Error: Could not open " << filename << endl;
        return 1;
    }

    vector<double> X, Y;
    double a, b;

    while (infile >> a >> b) {
        X.push_back(a);
        Y.push_back(b);
    }
    infile.close();

    int choice;
    cout << "\nChoose a task:\n";
    cout << "1. Print N lines\n";
    cout << "2. Calculate magnitudes\n";
    cout << "3. Least squares line fit\n";
    cout << "4. Calculate x^y (no loops/pow)\n";
    cout << "5. User magnitude (template)\n";
    cout << "Enter task number: ";
    cin >> choice;

    switch(choice) {
        case 1: task1(X, Y); break;
        case 2: task2(X, Y); break;
        case 3: task3(X, Y); break;
        case 4: task4(X, Y); break;
        case 5: task5(); break;
        default: cout << "Invalid choice.\n";
    }

    return 0;
}
