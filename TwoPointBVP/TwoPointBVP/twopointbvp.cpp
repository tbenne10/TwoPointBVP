#include "twopointbvp.h"
using namespace std;

TwoPointBVP::TwoPointBVP(double *dom, double (*dfunc)(vector<double> &)){
    domain = dom;
    diffusion = dfunc;
    reactionispresent = false;
    externalforceispresent = false;
    }

void TwoPointBVP::set_left_bdry(bool _leftisDirichlet, double *val){
    leftbdryisDirichlet = _leftisDirichlet;
    leftbdryvalues = val;
    }

void TwoPointBVP::set_right_bdry(bool _rightisDirichlet, double *val){
    rightbdryisDirichlet = _rightisDirichlet;
    rightbdryvalues = val;
     }

void TwoPointBVP::set_reaction(double(*func)(vector<double>&))
{
	reaction = func;
    reactionispresent = true;
}

void TwoPointBVP::set_reaction(double (*functone)(vector<double> &), double (*functtwo)(vector<double> &)){
    reaction = functone;
    partialreactionpartialu = functtwo;
    reactionispresent = true;
     }

void TwoPointBVP::set_externalforce(double (*func)(vector<double> &)){
    externalforce = func;
    externalforceispresent = true;
}

 double * TwoPointBVP::get_domain() const{
    return domain;
    }

 bool TwoPointBVP::left_bdry_is_Dirichlet() const{
    return leftbdryisDirichlet;
}

bool TwoPointBVP::right_bdry_is_Dirichlet() const{
    return rightbdryisDirichlet;
    }

double * TwoPointBVP::get_left_bdry_values() const{
    return leftbdryvalues;
    }

double * TwoPointBVP::get_right_bdry_values() const{
    return rightbdryvalues;
    }

 bool TwoPointBVP::reaction_is_present() const {
    return reactionispresent;
    }

bool TwoPointBVP::externalforce_is_present() const {
    return externalforceispresent;
}

 double TwoPointBVP::eval_diffusion(vector<double> &x) const{
     return diffusion(x);
    }

 vector<double> TwoPointBVP::eval_reaction(vector<double> &x) const{ 
	 double one, two;
	 one = reaction(x);
     two = partialreactionpartialu(x);
    vector<double> vec(2);
	vec[0] = one;
	vec[1] = two;
     return vec;
     // return reaction(x);
    }



 double TwoPointBVP::eval_reaction(vector<double> &x, int i) const {
	 return reaction(x);
 }

double TwoPointBVP::eval_externalforce(vector<double> &x) const{
    return externalforce(x);
}

TwoPointBVP::~TwoPointBVP(){
    ;
}
