#ifndef twopointbvp_h
#define twopointbvp_h


#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

// A C++ class definition for two−point BVP



class TwoPointBVP
 {
    protected:
    double * domain;
    bool leftbdryisDirichlet; //value is "true" if left boundary is Dirichlet
    bool rightbdryisDirichlet; //value is "true" if right boundary is Dirichlet
    double * leftbdryvalues; // contains gamma_0 and g_0 (see course notes).
    double * rightbdryvalues; // contains gamma_0 and g_0 (see course notes).
    bool reactionispresent; // value is "true" if there is reaction.
     bool externalforceispresent; // value is "true" if there is external force.
    double (*diffusion)(std::vector<double> &);// diffusion coefficient, k(x)
    double(*reaction)(std::vector<double> &);// reaction coefficient, r(x)
    double(*externalforce)(std::vector<double> &);//external force coefficient, f(x)
    double (*partialreactionpartialu) (vector<double> &); //Added from Q5
     
    public:
    
    // Constructor of the class, it takes two input
     TwoPointBVP(double *dom, double (*func)(vector<double> &));
    // Set data for the left boundary.
    void set_left_bdry(bool _leftisDirichlet, double *val);
    
    // Set data for the right boundary.
    void set_right_bdry(bool _rightisDirichlet, double *val);
    
    // Set reaction.
    void set_reaction(double(*func)(vector<double> &)); //OLD
    void set_reaction(double (*functone)(vector<double> &), double (*functtwo)(vector<double> &));
     //func 1 is for r(x, u) while func2 is for dr/du(x,u)
    
    //set ext. force
     void set_externalforce(double(*func)(vector<double> &));
     
    /***************
        Functions below can be called from outside TwoPointBVP to access various
        pertaining aspects of TwoPointBVP
     ***************/
    
    // Return domain (see above)
    double *get_domain() const;
    
    // Return leftbdryisDirichlet
    bool left_bdry_is_Dirichlet() const;
    
    // Return rightbdryisDirichlet
    bool right_bdry_is_Dirichlet() const;
    
    // Return leftbdryvalues, i.e., one that contains values of gamma_0 and g_0
    double * get_left_bdry_values() const;
    
    // Return leftbdryvalues, i.e., one that contains values of gamma_L and g_L
    double * get_right_bdry_values() const;
    
    // Return reactionispresent
    bool reaction_is_present() const;
     
     // return externalforcepresent;
    bool externalforce_is_present() const;
     
    // Return the value of k(x)
     double eval_diffusion(vector<double> &x) const;
     // Return the value of r(x)
     vector<double> eval_reaction(vector<double> &x) const;
	 double eval_reaction(vector<double> &x,int i) const;
     //return the value of f(x)
     double eval_externalforce(vector<double> &x) const;
     
     ~TwoPointBVP();
     };

#endif
