// -------------------------------------------------------------------
// -----       CosmicsGenerator source file for SHiP             -----
// -----       Version by 01.06.15  by Martin Franke             -----
// -----       mailto: mfranke(at)physik.hu-berlin.de            -----
// -------------------------------------------------------------------

#include "TROOT.h"
#include "FairPrimaryGenerator.h"
#include "CosmicsGenerator.h"
#include "TDatabasePDG.h"               // for TDatabasePDG
#include "TMath.h"

using namespace std;

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

Bool_t CosmicsGenerator::Init(Float_t zmiddle, Bool_t largeMom){
	//general
	fRandomEngine = new Co3Rng();
	TDatabasePDG* pdgBase = TDatabasePDG::Instance();
	mass = pdgBase->GetParticle(13)->Mass(); // muons!

	EVENTS = 400000; // #simulated events per "spill"
	// coordinate system
	xdist = 3000; // production area size [cm]
	zdist = 9000; // production area size [cm]
	yTop = 600; // box top layer [cm]-> also change weight3 accordingly when varying
	xTop = 300; // box side layer [cm]-> also change weight3 accordingly when varying
	zTop = 3650; // box length [cm]-> also change weight3 accordingly when varying
	z0 = zmiddle; // relative coordinate system [cm] (Z_muonstation + (Z_veto - 2 * Z_Tub1))/2,... Z_veto <0 ! ->z0 = 716
	cout<<"z0: "<<z0<<endl;
	
	high = largeMom;
	if (high) cout<<"Simulation for high momentum"<<endl;
	else cout<<"Simulation for low momentum"<<endl;
	
	// weights
	double weight1;
	if (!high) { // momentum range 1 GeV - 100 GeV
		weight1 = 123*xdist*zdist/EVENTS/10000; // expected #muons per spill/ #simulated events per spill: 123*30*90/500000
	} 
	else { // momentum range 100 GeV - 1000 GeV
		double I = fRandomEngine->fSpectrumH->Integral(100,1000);
		weight1 = 2*TMath::Pi()/3*I*xdist*zdist/EVENTS/10000; // expected #muons per spill/ #simulated events per spill: 123*30*90/500000
	}
	double weight3 = 4.834154338; // MC average of nTry/nEvents 4.834154338 +- 0.000079500												
	weight = weight1 / weight3;

	// running
	y = 1900; //all muons start 19m over beam axis
	nInside = 0;  nEvent = 0; nTest = 0; weighttest = 0; // book keeping
	return kTRUE;
}

// -----   Passing the event   ---------------------------------------------
Bool_t CosmicsGenerator::ReadEvent(FairPrimaryGenerator* cpg){
	Bool_t hit = 0;

	do{
		// shower characteristics
		double phi = fRandomEngine->Uniform(0,2*TMath::Pi());
		double theta = fRandomEngine->fTheta->GetRandom();

		//momentum components
		px = TMath::Sin(phi)*TMath::Sin(theta); 
		pz = TMath::Cos(phi)*TMath::Sin(theta);
		py = -TMath::Cos(theta);

		// start position, area 1120 m^2
		x = fRandomEngine->Uniform(-xdist/2,xdist/2);
		z = fRandomEngine->Uniform(z0 - zdist/2, z0 + zdist/2);

		// claim for flight through a box surrounding the detector
		if((abs(x-(y+yTop)*px/py) < xTop && abs(z-z0-(y+yTop)*pz/py) < zTop) || (abs(x-(y-yTop)*px/py) < xTop && abs(z-z0-(y-yTop)*pz/py) <  zTop)|| abs(y-(x+xTop)*py/px)<yTop && abs(z-z0-(x+xTop)*pz/px)<zTop || abs(y-(x-xTop)*py/px)<yTop && abs(z-z0-(x-xTop)*pz/px)<zTop){
			
			//muon or anti-muon
			if (fRandomEngine->Uniform(0,1) < 1.0/2.278){id = 13;}
			else{id = -13;}
			
			if (!high) P = fRandomEngine->fSpectrumL(theta);
			else P = fRandomEngine->fSpectrumH->GetRandom();
				
			px = px*P;
			py = py*P;
			pz = pz*P;
					
			// transfer to Geant4
			cpg->AddTrack(id,px,py,pz,x,y,z,-1,true,TMath::Sqrt(P*P+mass*mass),0,weight);  // -1 = Mother ID, true = tracking, SQRT(x) = Energy, 0 = t
			hit = 1; nInside++;
		}
		weighttest += weight;	nTest++;
	}while(!hit);
	nEvent++;
	if (!nEvent%10000){cout<<nEvent/10000<<"10k events have been simulated"<<endl;}
	return kTRUE;
}

// -------------------------------------------------------------------------

ClassImp(CosmicsGenerator)
