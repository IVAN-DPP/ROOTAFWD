#ifndef MAXLIKE_H
#define MAXLIKE_H
#include "TVirtualFitter.h"
#include "TStyle.h"
//#include "/usr/local/root/include/Minuit2/FCNBase.h"
//#include "TFitterMinuit.h"
#include "TSystem.h" 
#include "Math/Factory.h"
//#include "Minuit2/FCNBase.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/IFunction.h"
#include <vector>
#include <iostream>

using namespace std;

class MaxLike: public ROOT::Math::Functor{

private:
    vector<double> Phi;
    vector<double> GammaP;
public:
    MaxLike(vector<double>  INPhi,  vector<double> INGammaP){
        Phi=INPhi;
        GammaP=INGammaP;
	//std::cout << "Phi: " << Phi.at(0) << std::endl; 
	//std::cout << "Pol: " << GammaP.at(0) << std::endl; 
    }
    double operator() (const double * x) const {
        double sum=0;
        typedef vector<int>::size_type vec_sz;
        vec_sz m=Phi.size();
        for (vec_sz i=0; i<m; i++) {
            sum+=-TMath::Log(1-GammaP[i]*x[0]*TMath::Cos(2*Phi[i]+2*x[1]*TMath::DegToRad()));
	    
        }
        return sum;
    }
    double Up() const { return 0.5; }
};
#endif

