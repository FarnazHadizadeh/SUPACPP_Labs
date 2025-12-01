  C++ Assignment: 

This project contains a modular C++ program to read a data file of (x, y) pairs and perform multiple analysis tasks, including printing lines, calculating magnitudes, performing least squares line fitting, and computing powers without loops/pow.  

The program is designed to be **modular**: the **main program** is in `AnalyseData.cxx`, and all task-specific functions are implemented in `CustomFunctions.cxx` with declarations in `CustomFunctions.h`.

---

## Files

- `AnalyseData.cxx`  
  - Contains the **main function**.  
  - Reads the input data file `data.txt` into vectors.  
  - Displays a menu for the user to select which task to run.  
  - Calls the corresponding task function from `CustomFunctions.cxx`.  

- `CustomFunctions.cxx`  
  - Implements all task functions (Tasks 1–5):  
    1. Print first N lines of data.  
    2. Calculate magnitudes of vectors.  
    3. Perform least squares line fitting and chi-squared calculation.  
    4. Calculate x^y without using loops or `pow()`.  
    5. Calculate the magnitude of a user-entered vector using a template function.  

- `CustomFunctions.h`  
  - Header file declaring all task functions so they can be called from `AnalyseData.cxx`.  


