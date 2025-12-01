#include "CustomFunctions.h"
#include <iostream>
#include <fstream>
#include <cmath>

//testcomment

using namespace std;

// ------------------- Task 1 -------------------
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

// ------------------- Task 2 -------------------
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

// ------------------- Task 3 -------------------
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

// ------------------- Task 4 -------------------
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

// ------------------- Task 5 Template -------------------
template <typename T>
T magnitude(T x, T y) {
    return sqrt(x*x + y*y);
}

// Wrapper for Task 5
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

// Explicit instantiation of template for double
template double magnitude(double, double);