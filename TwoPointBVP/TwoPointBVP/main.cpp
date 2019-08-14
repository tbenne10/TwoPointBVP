//
//  main.cpp
//  
//
//  Created by Trystan Bennett on 5/5/18.
//  Copyright  2018 Trystan Bennett. All rights reserved.
//

#include <stdio.h>
#include <stdio.h>
#include "tridiagonal_matrix.h"
#include "twopointbvp.h"
#include "twopointbvpappr.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <complex>

using namespace std;



double diffusioncoeff(vector<double> &par){
    //P1
	return 1.0;

    //P2
    //return 1.0 / (1.0 - 0.99 * par[0] * sin(5.11 * 3.14159265 * par[0]));

    //P3
    //return 1.0;


}

double reactioncoeff(vector<double> &par) {
	//P1
	return -(100.0/88000000.0);

	//P2
	//return 0.0;

	//P3
	//double u = par[1];
	//double x = -pow(2.7182818, u);
	//return x;

}
double partialcoeff(vector<double> &par) {
	//P3
	double u = par[1];
	double x = -pow(2.7182818, u);
	return x;

}

double forcecoeff(vector<double> &par){
	//P1
	double x = -((200.0*par[0]) / (2.0*88000000.0));
	double y = (50.0 - par[0]);
	return x*y;

	//P2
	//return 1.0;

	//P3
	//return 0.0;

}

bool true_sol_is_present = true;

//Set this depending on if the reaction is semi linear
bool semiLinear = false; 

double truesol(double x){
	//P1
	double g = -1.0 / (1000.0*(1.0+exp(5.0*sqrt(5.0/11.0)/2.0)))	* 
		exp(-x/(4.0*sqrt(55.0)))*(-exp(x / (4.0*sqrt(55.0)))*(x*x - 50.0*x + 1760.0) -
			exp((x+50.0)/(4.0*sqrt(55.0)))*(x*x - 50.0*x + 1760.0)+ 1760.0*exp(x/(2.0*sqrt(55.0))) + 1760.0*exp((5.0*sqrt(5.0/11.0))/2.0));
	return g;


	//P2
	//return (x*x/2 + (-0.472297)*x - (0.000478578-0.0616686*x*x)*cos(16.0535 * x) - 
	//(0.00768287*x*sin(16.0535*x)) - 0.00384144*(-0.472297)*sin(16.053* x) +
		//0.0616686*(-0.472297)*x*cos(16.0535*x) + 0.000478578);

	//P3
	//double theta = 1.51716;
	//double theta = 10.9387;
	//return -2*log((cosh(.5*(x - .5)*theta) / cosh(.25*theta)));
}

int main()
{
    
    // Set the two point boundary value problem.
    //
    double *dom = new double[2]; //P1: domain 0 -> 50, P2: 0->1, P3: 0->1
    dom[0] = 0.0;
    dom[1] = 50.0;
    TwoPointBVP *prob = new TwoPointBVP(dom, diffusioncoeff);
    
    double *Lbval = new double[2]; 
    Lbval[0] = 0.0;//P1:0 P2:0 P3:0
    Lbval[1] = 0.0; //P1:0 P2:0 P3:0
    prob->set_left_bdry(true, Lbval);
    
    double *Rbval = new double[2];  
    Rbval[0] = 0.0; //P1:0 P2:0 P3:0
    Rbval[1] = 50.0; //P1:50 P2:1 P3:1
    prob->set_right_bdry(true, Rbval); 
    prob->set_externalforce(forcecoeff);

	if (semiLinear == true) {
		prob->set_reaction(reactioncoeff, partialcoeff);
	}
	else {
		prob->set_reaction(reactioncoeff);
	}
    
    // info regarding the BVP.
    //
    std::cout << "\n Several info regarding the two point BVP:\n\n";
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
    //
    int numsubintervals = 10;

	//Variables for Ui
	vector<double> sol;
    vector<double> xcoord;
    vector<double> er; //vector of errors to find the max of
    double * erH = new double[numsubintervals-1]; //Array of max error for each h
    double * sub = new double[numsubintervals-1]; //Array of H
    sub[0] = numsubintervals; //initial
    double * inRes = new double[numsubintervals];

	//Variables for U_hati 
	vector<double> sol2;
	vector<double> xcoord2;
	vector<double> erW; //vector of errors to find the max of
	double * erWH = new double[numsubintervals - 1]; //Array of max error for each h
	double * sub2 = new double[(2 * numsubintervals) - 1]; //Array of H
	sub2[0] = 2*numsubintervals; //initial
    
    
    for (int i = 1; i < 10; i++) {
        sub[i] = sub[i - 1] * 2;
		sub2[i] = sub[i] * 2;
    }

    for (int i = 0; i < 10; i++) {
        numsubintervals = sub[i];
        double * subintervals = new double[numsubintervals];
		double * subintervals2 = new double[numsubintervals*2];
        for (int i = 0; i < numsubintervals; i++)
            subintervals[i] = (dom[1] - dom[0]) / numsubintervals;
        inRes[i] = subintervals[i];
        
		for (int i = 0; i < numsubintervals*2; i++) {
			subintervals2[i] = (dom[1] - dom[0]) /(2*numsubintervals);
		}
        TwoPointBVPAppr * method = new TwoPointBVPAppr(numsubintervals, subintervals, prob);
		TwoPointBVPAppr * method2 = new TwoPointBVPAppr(numsubintervals*2, subintervals2, prob);
		if (semiLinear == true) {
			sol = method->Solve(10, .1);
			sol2 = method2->Solve(10, .1);

		}
		else {
			sol = method->Solve();
			sol2 = method2->Solve();
		}

		//Loop here to set Wi
		vector<double> W;
		for (int i = 0; i < numsubintervals+1; i++) {
			W.push_back(.33*((4 * sol2[2 * i]) - sol[i]));
		}
        xcoord = method->get_xcoord();
		xcoord2 = method2->get_xcoord();
        if (true_sol_is_present)
        {
            double s = (dom[1] - dom[0]) / numsubintervals;
			double s2 = (dom[1] - dom[0]) / (numsubintervals*2);
            double x = dom[0];
            for (int i = 0; i < numsubintervals + 1; i++)
            {
				//cout << truesol(x) << "    " << sol.at(i) << "    " << W[i] <<  endl;
                er.push_back(abs(sol.at(i) - truesol(x))); //err =approx slt - true sol
				erW.push_back(abs(W[i] - truesol(x)));
                x += s;
            }
            //fileout.close();
            auto error = std::max_element(std::begin(er), std::end(er));
            double max = *error;
			erH[i] = max;
			auto error2 = std::max_element(std::begin(erW), std::end(erW));
			max = *error2;
			erWH[i] = max;
			er.clear();
			erW.clear();
        }
        delete method; //delete to refresh
        delete[]subintervals;
    }
    //end loop
    ofstream fileout;
    fileout.open("p3_2.txt");
    for (int i=0; i<10; i++){
        fileout << inRes[i] << " " << erWH[i] << " " << endl;
    }
    fileout.close();
	std::cout << "The following error values for Ui, Wi at each stepsize" << endl;
	for (int i = 0; i < 10; i++) {
		std::cout << "Step size: "  <<inRes[i] << std::setw(20) << "Ui: " << erH[i] << std::setw(20) << "Wi:  " << erWH[i] <<endl;
	}


    delete []Rbval;
    delete []Lbval;
    delete prob;
    delete []dom;
    system("pause");
    return 0;
}
