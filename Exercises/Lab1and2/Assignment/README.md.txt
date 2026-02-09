# SUPAC++ Assignment 2 – Labs 3 & 4

**Author:** Farnaz  
**Course:** SUPAC++  
**Assignment:** Statistical Distributions and Sampling  
**Due date:** 12th December

---

## Overview

This assignment investigates the probability distribution used to generate an unknown one-dimensional dataset.
A provided executable generates mystery data sampled from an unknown statistical distribution.

The aim is to:
- Compare mystery data with analytical distributions
- Implement these distributions using C++ classes and inheritance
- Numerically integrate the functions
- Identify the most likely underlying distribution
- (Optionally) sample new data using the Metropolis algorithm

---

## Project Structure
SUPAC_Assignment2/
│
├── src/
│ ├── FiniteFunctions.h
│ ├── FiniteFunctions.cxx
│ ├── CustomDistributions.h
│ ├── CustomDistributions.cxx
│ └── main.cxx
│
├── Output/
│ └── Data/
│ └── MysteryDataXXXXX.txt
│
├── Makefile
└── README.md
