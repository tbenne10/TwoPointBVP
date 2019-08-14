//
//  twopointbvpappr.h
//  BVPsolve//
//  Created by Trystan Bennett on 3/3/18.
//  Copyright  2018 Trystan Bennett. All rights reserved.
//

#ifndef twopointbvpappr_h
#define twopointbvpappr_h

#include <stdio.h>

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "twopointbvp.h"
#include "tridiagonal_matrix.h"
using namespace std;

class TwoPointBVPAppr
{
    protected:
    int numsubintervals;
    const double *steplengths;
    const TwoPointBVP * theproblem;
    vector<double> xcoord;
    vector<double> midcoord;
    vector<double> Deltax;
    
    void AssembleDiffusion(tridiagonal_matrix * tmat);
	vector<double> AssembleReaction();
    void AssembleReaction(vector<double> &a, vector<double> &b, vector<double> &c);
    
    //vector<double> AssembleReaction();
    
    vector<double> AssembleForce();
    
    
    public:
    
    TwoPointBVPAppr(int N, const double * subintervallengths, const TwoPointBVP * prob);
    
    vector<double> get_xcoord();
    
    vector<double> Solve(int max_num_iter, double TOL);
	vector<double> Solve();
    ~TwoPointBVPAppr();
    };



#endif /* twopointbvpappr_hpp */
