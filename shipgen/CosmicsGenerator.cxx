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
	if (!high) { // momentum range 1 GeV - 100 GeV
		weight1 = 123*xdist*zdist/EVENTS/10000; // expected #muons per spill/ #simulated events per spill: 123*30*90/500000
		weight2 = TMath::Pi/3;
	} 
	else {
		double I = fRandomEngine->fspectrumH->Integral(100,1000);
		weight1 = 2*TMath::Pi()/3*I*xdist*zdist/EVENTS/10000; // expected #muons per spill/ #simulated events per spill: 123*30*90/500000
		weight2 = 900/I; // 1/(mean momentum weight), P_max-P_min/(3*0.3044/2pi)
	}
	double weight3 = 4.834154338; // MC average of nTry/nEvents 4.834154338 +- 0.000079500												
	weight = weight1 * weight2  / weight3;

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
			
			if (!high) { P = fRandomEngine->fSpectrumL(theta); w = weight;}
			else { 
				P = fRandomEngine->Uniform(100,1000);
				w = weight * fRandomEngine->fspectrumH->Eval(P);
			}
			px = px*P;
			py = py*P;
			pz = pz*P;
					
			// transfer to Geant4
			cpg->AddTrack(id,px,py,pz,x,y,z,-1,true,TMath::Sqrt(P*P+mass*mass),0,w);  // -1 = Mother ID, true = tracking, SQRT(x) = Energy, 0 = t
			hit = 1; nInside++;
		}
		weighttest += w;	nTest++;
	}while(!hit);
	nEvent++;
	if (!nEvent%10000){cout<<nEvent/10000<<"10k events have been simulated"<<endl;}
	return kTRUE;
}

// -------------------------------------------------------------------------

ClassImp(CosmicsGenerator)
