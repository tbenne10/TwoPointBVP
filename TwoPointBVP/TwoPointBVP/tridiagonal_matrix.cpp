//
//  tridiagonal_matrix.cpp
//  TriDiag
//
//  Created by Trystan Bennett on 2/21/18.
//  Copyright  2018 Trystan Bennett. All rights reserved.
//

#include "tridiagonal_matrix.h"

tridiagonal_matrix::tridiagonal_matrix(int m)
{
    dimension = m;
    diag.resize(m);
    upperdiag.resize(m-1);
    lowerdiag.resize(m-1);
	r.resize(m-1);
	hatupperdiag.resize(m-1);
    transformed = false;
    }

tridiagonal_matrix::tridiagonal_matrix(tridiagonal_matrix *mat)
{
	dimension = mat->get_dimension();
	diag.resize(dimension);
	upperdiag.resize(dimension - 1);
	lowerdiag.resize(dimension - 1);
	r.resize(dimension - 1);
	hatupperdiag.resize(dimension - 1);
	for (int i = 0; i<dimension; i++)
		diag[i] = mat->get_diagonal_entry(i);
	for (int i = 0; i<dimension - 1; i++)
	{
		upperdiag[i] = mat->get_upper_diagonal_entry(i);
		lowerdiag[i] = mat->get_lower_diagonal_entry(i);
	}
	transformed = mat->istransformed();
	if (transformed)
	{
		for (int i = 0; i<dimension - 1; i++)
		{
			r.push_back(mat->get_r_entry(i));
			hatupperdiag.push_back(mat->get_hat_upper_diagonal_entry(i));
		}
		transformed = true;
	}
}


int tridiagonal_matrix::get_dimension()
{
    return dimension;
    }

void tridiagonal_matrix::set_diagonal_entry(int i, double val) {
    diag[i] = val;
    }

void tridiagonal_matrix::set_upper_diagonal_entry(int i, double val)
{
    upperdiag[i] = val;
    }
void tridiagonal_matrix::set_lower_diagonal_entry(int i, double val) {
    lowerdiag[i] = val;
    }

double tridiagonal_matrix::get_diagonal_entry(int i)
{
    return diag[i];
    }

double tridiagonal_matrix::get_upper_diagonal_entry(int i){
    return upperdiag[i];
    }
double tridiagonal_matrix::get_lower_diagonal_entry(int i)
{
    return lowerdiag[i];
    }

double tridiagonal_matrix::get_hat_upper_diagonal_entry(int i){
    return hatupperdiag[i];
}

double tridiagonal_matrix::get_r_entry(int i){
    return r[i];
}

bool tridiagonal_matrix::istransformed(){
    return transformed;
    }

void tridiagonal_matrix::add_to_diagonal_entry(int i, double val)
{
    diag[i] += val;
    }

void tridiagonal_matrix::add_to_upper_diagonal_entry(int i, double val)
{
    upperdiag[i] += val;
    }

void tridiagonal_matrix::add_to_lower_diagonal_entry(int i, double val){
    lowerdiag[i] += val;
}

// This is where the matrix is transformed into upper triangular.
void tridiagonal_matrix::transform()
{
   //  Implement here the transformation. Fill the entries of r and hatupperdiag.

	//r is f and e is hatupperdiag
	
	hatupperdiag[0] = upperdiag[0] / diag[0];
	
	for (int i = 1; i <= dimension - 2; i++) {
		r[i - 1] = diag[i] - (lowerdiag[i - 1] * hatupperdiag[i - 1]);
		hatupperdiag[i] = upperdiag[i] / r[i - 1];
	}
	r[dimension - 2] = diag[dimension - 1] - ((lowerdiag[dimension - 2]*hatupperdiag[dimension-2]));

	 transformed = true;
    }

// Given vector rhs, return sol satisfying A sol = rhs.
vector <double> tridiagonal_matrix::solve_linear_system(const vector<double> & rhs)
{
    vector<double> hatrhs(dimension); // modified right hand side
    vector<double> sol(dimension); // solution
	transform();
	hatrhs[0] = rhs[0] / diag[0];
	for (int i = 1; i <= dimension - 2; i++) {
		hatrhs[i] = (rhs[i] - (lowerdiag[i-1]*hatrhs[i-1]))/r[i-1];
	}
	hatrhs[dimension-1] = ((rhs[dimension-1] - (lowerdiag[dimension - 2] * hatrhs[dimension - 2]))) / r[dimension - 2];


	//backward solve
	sol[dimension - 1] = hatrhs[dimension - 1];
	for (int j = dimension - 2; j >= 0; j--) {
		sol[j] = hatrhs[j] - (hatupperdiag[j]*sol[j+1]);
	}


    return sol;
}

// Given vector lhs, return rhs = A lhs.
vector <double> tridiagonal_matrix::Mult(const vector<double > & lhs)
{
	vector <double > rhs(dimension, 0.0);
	//rhs[0] = (diag[0] + upperdiag[0])*lhs[0];
	rhs[0] = diag[0] * lhs[0] + upperdiag[0] * lhs[1];
	for (int i = 1; i <= dimension - 2; i++) {
		//rhs[i] = (lowerdiag[i] + diag[i] + upperdiag[i])*lhs[i];
		rhs[i] = lowerdiag[i - 1] * lhs[i - 1] + diag[i] * lhs[i] + upperdiag[i] * lhs[i + 1];
	}
	//rhs[dimension-1] = (lowerdiag[dimension-2] + diag[dimension-1])*lhs[dimension-1];
	rhs[dimension - 1] = lowerdiag[dimension - 2] * lhs[dimension - 2] + diag[dimension - 1] * lhs[dimension - 1];
	return rhs;
    }
// Empty destructor
tridiagonal_matrix::~tridiagonal_matrix(){
    ;}
