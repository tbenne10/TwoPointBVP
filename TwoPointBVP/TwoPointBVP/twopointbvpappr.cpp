//
//  twopointbvpappr.cpp
//  BVPsolver
//
//  Created by Trystan Bennett on 3/3/18.
//  Copyright  2018 Trystan Bennett. All rights reserved.
//


#include "twopointbvpappr.h"
#include "tridiagonal_matrix.h"
#include "twopointbvp.h"
#include "twopointbvpappr.h"
#include <algorithm>


TwoPointBVPAppr::TwoPointBVPAppr(int N, const double * subintervallengths, const TwoPointBVP * prob){
    numsubintervals = N;
    steplengths = subintervallengths; //h_i
    theproblem = prob;
    // If either reaction or externalforce is present, then we need Delta x_i
    if (theproblem->reaction_is_present() || theproblem->externalforce_is_present()){
            Deltax.resize(numsubintervals+1);
            Deltax[0] = 0.5 * steplengths[0];
            for (int i=1; i<numsubintervals; i++)
                Deltax[i] = 0.5 * (steplengths[i-1]+steplengths[i]);
                Deltax[numsubintervals] = 0.5 * steplengths[numsubintervals-1];
               }
 
   // This is the x_i
    xcoord.resize(numsubintervals+1);
    double * domain = theproblem->get_domain();
    xcoord[0] = domain[0];
    for (int i=1; i<=numsubintervals; i++)
        xcoord[i] = xcoord[i-1] + steplengths[i-1];
        midcoord.resize(numsubintervals+2);
        midcoord[0] = domain[0];
    for (int i=1; i<=numsubintervals; i++){
            midcoord[i] = 0.5 * (xcoord[i] + xcoord[i-1]);
        
    }
    midcoord[N+1] = xcoord[numsubintervals];
            }
 vector<double> TwoPointBVPAppr::get_xcoord() {
   return xcoord;
    }



//** DIFUSION **
//Fills out the tridiagonal matrix
//Diffusion coefficient is k

void TwoPointBVPAppr::AssembleDiffusion(tridiagonal_matrix * tmat)
{
    vector <double> kappa(numsubintervals+2);
    vector<double> par(1);
    for (int i=0; i<numsubintervals+2; i++)
        {
            par[0] = midcoord[i];
            kappa[i] = theproblem->eval_diffusion(par);
            }
    
    for (int i=0; i< numsubintervals + 2; i++)
        kappa[i] = kappa[i] / steplengths[i-1];
    
        // taking care of zeroth row, i.e., left boundary.
        if (theproblem->left_bdry_is_Dirichlet())
            {
                
                 tmat->set_diagonal_entry(0, 1);
                 tmat->set_upper_diagonal_entry(0, 0);
                }
    else
        {
            tmat->set_upper_diagonal_entry(0, -kappa[1]);
            tmat->set_diagonal_entry(0, kappa[1] - theproblem->get_left_bdry_values()[0]);
            
            }
    // filling up the matrix tmat for internal row.
    for (int i=1; i<numsubintervals; i++)
        {
            tmat->set_diagonal_entry(i, (kappa[i]+kappa[i+1]));
            tmat->set_upper_diagonal_entry(i, -kappa[i+1]);
            tmat->set_lower_diagonal_entry(i-1, -kappa[i]);
            
             }
   
   // taking care of the last row, i.e., right boundary.
     if (theproblem->right_bdry_is_Dirichlet())
        {
            tmat->set_diagonal_entry(numsubintervals, 1);
            tmat->set_lower_diagonal_entry(numsubintervals-1, 0);
            }
    else
        {
            
            tmat->set_lower_diagonal_entry(numsubintervals-1, -kappa[numsubintervals]);
            tmat->set_diagonal_entry(numsubintervals, kappa[numsubintervals] - theproblem->get_right_bdry_values()[0]);
            
        }
   }

vector<double> TwoPointBVPAppr::AssembleReaction()
{
	// This vector contains the reaction algebraic terms
		 vector<double> Reac(numsubintervals + 1);

		// taking care of zeroth entry, i.e., left boundary.
		 if (theproblem->left_bdry_is_Dirichlet())
		{
		 // Add zero to zeroth row
			Reac[0] = 0;
		 }
	 else
		{
		// Add f(x_0) * deltax_0
		 vector<double> x_reaction(1);
		 x_reaction[0] = xcoord[0];
		
			 Reac[0] = theproblem->eval_reaction(x_reaction,1) * Deltax[0];
		 }
	
		 for (int i = 1; i < numsubintervals; i++)
		 {
		 // Add f(x_i) * deltax_i
			 vector<double> x_reaction(1);
		 x_reaction[0] = xcoord[i];
		
			 Reac[i] = theproblem->eval_reaction(x_reaction,1) * Deltax[i];
		 }

		 // taking care of last entry, i.e., right boundary.
		 if (theproblem->right_bdry_is_Dirichlet())
		 {
		 // Add zero to last entry
			Reac[numsubintervals] = 0;
		 }
        else
		{
		 // Add f(xi) * deltaxi
		 vector<double> x_reaction(1);
		x_reaction[0] = xcoord[numsubintervals];
		
		 Reac[numsubintervals] = theproblem->eval_reaction(x_reaction,1) * Deltax[numsubintervals];
		 }
	
		 return Reac;
	 }

// ** REACTION **
//This is m+1 size, since we include the boundary conditions.
//No contribution from reaction if it is Dirichlet
void TwoPointBVPAppr::AssembleReaction(vector<double> &U, vector<double> &R, vector<double> &Rp)
{
    
    vector<double> Reac;
	vector<double> par(1);
	vector<double> xu_reaction(2);
    //taking care of zeroth entry, i.e., left boundary.
    if (theproblem->left_bdry_is_Dirichlet())
    {
    R[0] = 0.0;
    Rp[0] = 0.0;
    }
    else
    {
	
	xu_reaction[0] = xcoord[0];
	xu_reaction[1] = U[0];
	Reac = theproblem->eval_reaction(xu_reaction);
    R[0] = Deltax[0] * Reac[0];
    Rp[0] =Deltax[0] * Reac[1];
    }
    
    if (theproblem->right_bdry_is_Dirichlet())
    {
        R[numsubintervals] = 0.0;
        Rp[numsubintervals] = 0.0;
    }
    else
    {
		xu_reaction[0] = xcoord[numsubintervals];
		xu_reaction[1] = U[numsubintervals];
		Reac = theproblem->eval_reaction(xu_reaction);
        R[numsubintervals] = Deltax[numsubintervals] * Reac[0];
        Rp[numsubintervals] =Deltax[numsubintervals] * Reac[1];
    }
    
    for (int i=1; i<numsubintervals; i++){
		xu_reaction[0] = xcoord[i];
		xu_reaction[1] = U[i];
		Reac = theproblem->eval_reaction(xu_reaction);
        R[i] = Deltax[i] * Reac[0];
        Rp[i] =Deltax[i] * Reac[1];
        
    }
	return;
}

// **FORCE**
//This is m+1 size, since we include the boundary conditions.
// **FORCE**
//This is m+1 size, since we include the boundary conditions.
vector<double> TwoPointBVPAppr::AssembleForce() {
	// This vector contains the force algebraic terms.
	vector<double> FF(numsubintervals + 1);
	vector<double> par(1);
	// taking care of zeroth entry, i.e., left boundary.
	if (theproblem->left_bdry_is_Dirichlet())
	{

		FF[0] = theproblem->get_left_bdry_values()[1];

	}
	else
	{
		if (theproblem->externalforce_is_present()) {
			par[0] = xcoord[0];

			FF[0] = Deltax[0] * theproblem->eval_externalforce(par) - theproblem->get_left_bdry_values()[1];
		}
		else {
			FF[0] = -1 * theproblem->get_left_bdry_values()[1];
		}
	}
	if (theproblem->externalforce_is_present()) {
	for (int i = 1; i < numsubintervals; i++)
		{
		par[0] = xcoord[i];

		FF[i] = Deltax[i] * theproblem->eval_externalforce(par);
		}
	}
	// taking care of last entry, i.e., right boundary.
	if (theproblem->right_bdry_is_Dirichlet())
	{
		FF[numsubintervals] = theproblem->get_left_bdry_values()[1];
	}
	else
	{
		if (theproblem->externalforce_is_present()) {
			par[0] = xcoord[numsubintervals];
			FF[numsubintervals] = Deltax[numsubintervals] * theproblem->eval_externalforce(par) - theproblem->get_right_bdry_values()[1];
		}
		else {
			FF[numsubintervals] = -1.0 * theproblem->get_right_bdry_values()[1];
		}
		}

	return FF;
}


vector<double> TwoPointBVPAppr::Solve(int max_num_iter, double TOL)
{
    tridiagonal_matrix * Gp, * A;
	int s = 0; 
    //vector<double> R, Rp, F;
	vector<double> R(numsubintervals + 1, 0.0);
	vector<double> Rp(numsubintervals + 1, 0.0);
	vector<double> F(numsubintervals + 1, 0.0);
    A = new tridiagonal_matrix(numsubintervals+1);
    AssembleDiffusion(A);
    F = AssembleForce();
    vector<double> U(numsubintervals+1,2.0);
    vector<double> h(numsubintervals+1,0.0);
    vector<double> G(numsubintervals+1,0.0);
    vector<double> TestConv(numsubintervals+1,0.0);
    U[0] = F[0];
    double norm = 0;
    U[numsubintervals] = F[numsubintervals];
    // The iteration
	for (int iter = 1; iter <= max_num_iter; iter++)
	{
		// Copy A to Gp
		Gp = new tridiagonal_matrix(A);
		if (theproblem->reaction_is_present())
		{
			AssembleReaction(U, R, Rp);
			for (int i = 0; i < numsubintervals; i++)
				Gp->add_to_diagonal_entry(i, Rp[i]);
		}
		// Complete the rest of the code below to:
		// fill G to contain AU + R  F,
		vector<double> AU(numsubintervals , 0.0);
		AU = A->Mult(U);
		for (int i = 0; i <= numsubintervals; i++) {
			G[i] = AU[i] + R[i] - F[i];
		}

		// - solve the linear system Gp h = G,
		Gp->transform();
		h = Gp->solve_linear_system(G);
		// - update U based on h,
      

		for (int i = 0; i < numsubintervals; i++) {
			U[i] = U[i] - h[i];
		}

		for (int i = 0; i < numsubintervals; i++) {
			h[i] = fabs(h[i]);
		}
		norm = *std::max_element(h.begin(), h.end());

		if (norm < TOL) {
			if (iter == max_num_iter) {
				std::cout << "Max iterations reached. " << endl;
			}
			s = iter;
			delete Gp;
			break;
		}
	}
        // has been exceeded, hence convergence fails.
        //delete Gp;
   
    delete A;
	std::cout << "Iterations to convergence: " << s << endl;
    return U;
}

vector<double> TwoPointBVPAppr::Solve()
{
    tridiagonal_matrix *A, *LMat;
	vector<double> RJ, F;
	A = new tridiagonal_matrix(numsubintervals + 1);
	
	// Calculate the tridiagonal matrix coming from diffusion component.
	 AssembleDiffusion(A);
	
		 // Add reaction term to the diagonal entries of A if it is present.
		if (theproblem->reaction_is_present())
		{
		RJ = AssembleReaction();
		for (int i = 0; i < numsubintervals + 1; i++)
		 {
		   A->add_to_diagonal_entry(i, RJ[i]);
		 }
		}

	   // The first and last entries of F take care of boundary values.
		if (theproblem->externalforce_is_present())
		{
		F = AssembleForce();
		}
	
		 // Solve the system and return it.
	  return A->solve_linear_system(F);
	 }

 TwoPointBVPAppr::~TwoPointBVPAppr() {
    ;
    }
