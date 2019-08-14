//
//  testappr.cpp
//  BVPsolver
//
//  Created by Trystan Bennett on 3/3/18.
//  Copyright © 2018 Trystan Bennett. All rights reserved.
//

#include <stdio.h>
#include "tridiagonal_matrix.h"
#include "twopointbvp.h"
#include "twopointbvpappr.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <complex>

using namespace std;
//double lamb = .1;
//double theta = .450047;
//double theta = 17.4217;
//double lamb = .5;
//double theta = 1.03557;
//double theta = 13.0832;
//double lamb = 1;
//double theta = 1.51716;
//double theta = 10.9387;
//double lamb = 2;
//double theta = 2.35755;
//double theta = 8.5072;

double diffusioncoeff(vector<double> &par){
    //p1
	//return 1.0/(1.0+par[0]*par[0]);
	//p2
	return 1.0;
}

double reactioncoeff(vector<double> &par) {
	double x = 6.0*pow(par[1], .333333333333333);
	return x;
	//return -lamb;
}
double partialcoeff(vector<double> &par) {
	double u = par[1];
	double x = 2.0*pow(u, -.6666666666666666667);
	return x;
	//return -lamb;
}

double forcecoeff(vector<double> &par){
    return -12.0*par[0]*par[0];
	//return 1;
    }

bool true_sol_is_present = true;


double truesol(double x){
    return pow((1+(x*x)),3);
	//return -2*log((cosh(.5*(x - .5)*theta) / cosh(.25*theta)));
    }

int main()
{

    // Set the two point boundary value problem.
    //−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    double *dom = new double[2];
    dom[0] = 0.0;
    dom[1] = 1.0;
    TwoPointBVP *prob = new TwoPointBVP(dom, diffusioncoeff);

    double *Lbval = new double[2];
    Lbval[0] = 0.0;//P1: 0 P2: 0
    Lbval[1] = 1.0; //P1: 1 P2: 0
    prob->set_left_bdry(true, Lbval); //P1 TRUE P2 TRUE

    double *Rbval = new double[2];
    Rbval[0] = 0.0; //P1: 0 P2:0
    Rbval[1] = -12.0; //P1: -12 P2 0
    prob->set_right_bdry(false, Rbval); //P1 FALSE P2 TRUE
    prob->set_externalforce(forcecoeff);
	prob->set_reaction(reactioncoeff,partialcoeff);

    // Info regarding the BVP.
    //−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    std::cout << "\n Info regarding the two point BVP:\n\n";
    double *domain = prob->get_domain();
    cout << " Domain is (" << domain[0] << "," << domain[1] << ")" << endl;

    double *val;
    val = prob->get_left_bdry_values();
    if (prob->left_bdry_is_Dirichlet())
    cout << " Left boundary is Dirichlet with value " << val[1] << endl;
    else
    cout << " Left boundary is NOT Dirichlet with gamma = " << val[0] << " and g = "
    << val[1] << endl;

    val = prob->get_right_bdry_values();
    if (prob->right_bdry_is_Dirichlet())
    cout << " Right boundary is Dirichlet with value " << val[1] << endl;
    else
    cout << " Right boundary is NOT Dirichlet with gamma = " << val[0] << " and g = "
    << val[1] << endl;

    if (prob->reaction_is_present())
    {
        cout << " Reaction is present\n";
        }
    else
    cout << " No reaction\n";

    if (prob->externalforce_is_present())
    {
        cout << " External force is present\n";
        }
    else
    cout << " No external force\n";

    // Run the approximation.
    //−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    int numsubintervals = 10;
	vector<double> sol;
	vector<double> xcoord;
	vector<double> er; //vector of errors to find the max of
	double * erH = new double[9]; //Array of max error for each h
	double * sub = new double[9]; //Array of H
	sub[0] = 10; //initial
    double * inRes = new double[numsubintervals];


	for (int i = 1; i < 10; i++) {
		sub[i] = sub[i - 1] * 2;
	}
		//loop here
	for (int i = 0; i < 10; i++) {
		numsubintervals = sub[i];
		double * subintervals = new double[numsubintervals];
		for (int i = 0; i < numsubintervals; i++)
			subintervals[i] = (dom[1] - dom[0]) / numsubintervals;
            inRes[i] = subintervals[i];

		TwoPointBVPAppr * method = new TwoPointBVPAppr(numsubintervals, subintervals, prob);
		sol = method->Solve(30, .001);
		xcoord = method->get_xcoord();

		if (true_sol_is_present)
		{
			double s = (dom[1] - dom[0]) / numsubintervals;
			double x = dom[0];
			for (int i = 0; i < numsubintervals + 1; i++)
			{
				//std::cout << truesol(x) <<endl;
				er.push_back(sol.at(i) - truesol(x)); //err =approx slt - true sol
				x += s;
			}
			//fileout.close();
			auto error = std::max_element(std::begin(er), std::end(er));
			double max = *error;
			erH[i] = max;

		}
		delete method; //delete to refresh
		delete[]subintervals;
	}
	//end loop
	ofstream fileout;
	fileout.open("test.txt");
    for (int i=0; i<10; i++){
     fileout << inRes[i] << " " << erH[i] << " " << endl;
    }
	
	fileout.close();

    delete []Rbval;
    delete []Lbval;
    delete prob;
    delete []dom;
	//system("pause");
    return 0;
    }
