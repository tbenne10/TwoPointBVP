//
//  tridiagonal_matrix.h
//  TriDiag
//
//  Created by Trystan Bennett on 2/21/18.
//  Copyright  2018 Trystan Bennett. All rights reserved.
//

#ifndef tridiagonal_matrix_h
#define tridiagonal_matrix_h

#include <stdio.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

class tridiagonal_matrix {
    private:
    int dimension;
    vector <double> diag;
    vector <double> upperdiag;
    vector <double> lowerdiag;
    vector <double> hatupperdiag; // modified upper diagonal entries
    vector <double> r;
    bool transformed;
    public:
    tridiagonal_matrix(int m);
    tridiagonal_matrix(tridiagonal_matrix *mat);
    int get_dimension();
    void set_diagonal_entry(int i, double val);
    void set_upper_diagonal_entry(int i, double val);
    void set_lower_diagonal_entry(int i, double val);
    double get_diagonal_entry(int i);
    double get_upper_diagonal_entry(int i);
    double get_lower_diagonal_entry(int i);
    bool istransformed();
    // diag[i] = diag[i] + val;
    void add_to_diagonal_entry(int i, double val);
    // upperdiag[i] = upperdiag[i] + val;
    void add_to_upper_diagonal_entry(int i, double val);
    // lowerdiag[i] = lowerdiag[i] + val;
    void add_to_lower_diagonal_entry(int i, double val);
    
    //functions added in this assignment
    double get_r_entry(int i);
    double get_hat_upper_diagonal_entry(int i);
    
    vector <double> solve_linear_system(const vector<double> & rhs);
    
    void transform();// perform matrix vector multiplication
    vector<double > Mult(const vector<double > & lhs);
    ~tridiagonal_matrix();
    };

#endif /* tridiagonal_matrix_h */
