// -------------------------------------------------------------------
// -----       CosmicsGenerator source file for SHiP             -----
// -----       Version by 01.06.15  by Martin Franke             -----
// -----       mailto: mfranke(at)physik.hu-berlin.de            -----
// -------------------------------------------------------------------

#ifndef PNDCoGENERATOR_H
#define PNDCoGENERATOR_H 1

#include "TROOT.h"
#include "FairGenerator.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"
#include "TH1.h"

using namespace std;
class FairPrimaryGenerator;

class Co3Rng{
	public:
	   Co3Rng() {
			rng = new TRandom3(gRandom->GetSeed());
			fTheta = new TF1("f2","cos(x)*cos(x)",0,1.3734); // theta_zenith up to 78.7deg
			fTheta->SetNpx(10);
			fSpectrumH = new TF1("f4","1400*TMath::Power(x,-2.7)*(1/(1+x/115)+0.054/(1+x/850))",100,1000); // momentum above 100GeV
		};
	   virtual ~Co3Rng() {delete rng; delete fTheta; delete fSpectrumH;};
	   double Uniform(Float_t min, Float_t max){return rng->Uniform(min,max);};
	   TF1 *fSpectrumH;
	   TF1 *fTheta;
	   double fSpectrumL(double theta); // momentum below 100GeV
	private:
	   TRandom3 *rng; //!
};

double Co3Rng::fSpectrumL(double theta){
	//returns a random P between 1GeV and 100GeV taken from the 
	// zenith angle dependend momentum distribution from 
	// doi: 10.1016/j.nuclphysbps.2005.07.056.
	// Here, the inverse of the function is computed and a random number 
	// between 0 and 1 mapped to the interval [1, 100[ GeV
	
	theta = 180*theta/TMath::Pi(); // theta in degrees
	double a = -0.8816/10000 /(1/theta  -0.1117/1000 * theta) - 0.1096 - 0.01966*TMath::Exp(-0.02040*theta);
	double b = 0.4169/100 /(1/theta  -0.9891/10000 * theta) + 4.0395 - 4.3118*TMath::Exp(0.9235/1000*theta);

	double gamma = sqrt(-TMath::Ln10()*a);
	double offset = 0.5*(b  + 1/TMath::Ln10())/a;
	double norm = TMath::Erf(gamma*(TMath::Log(100)+offset)) - TMath::Erf(gamma*offset);
	
	double r3 = rng->Uniform();
	
	return exp(TMath::ErfInverse(r3*norm+TMath::Erf(gamma*offset))/gamma-offset);
}

class CosmicsGenerator : public FairGenerator{
 public:
  	/** constructor,destructor **/
	CosmicsGenerator(){};  
	virtual ~CosmicsGenerator(){
		delete fRandomEngine; 
		cout<<nEvent<<" events have been generated."<<endl;
		cout<<"There is a total of "<<nInside<<"/"<<nTest<<" muons that passed close enough to the detector."<<endl;
		cout<<"Including the given weight this corresponds to ";
		if (high) cout<<weighttest/xdist/zdist*10000/0.3044<<" spills (1 spill = "<<xdist*zdist*0.3044/10000;
		else cout<<weighttest/xdist/zdist*10000/123<<" spills (1 spill = "<<xdist*zdist*123/10000;
		cout<<" real cosmic muons = "<<EVENTS<<" simulated events)."<<endl;
	};
  
	/** public method ReadEvent **/
	Bool_t ReadEvent(FairPrimaryGenerator*);  //!
	//  virtual Bool_t Init(); //!
	virtual Bool_t Init(Float_t zmiddle, Bool_t largeMom); //!
  
 private:
	Co3Rng *fRandomEngine;//!
  
 protected:
	double px,py,pz,x,y,z,w, weighttest, weight, z0, mass, yTop,xTop,zTop,xdist, zdist;
	int id,nInside,nEvent,nTest,EVENTS;//!
	Bool_t high;
		
	ClassDef(CosmicsGenerator,3);
};

#endif /* !PNDCoGENERATOR_H */
