#include <iostream>
#include <cmath>
#include "Gobbiarray.h"
#include "frag.h"
#include "decay.h"
#include "constants.h"
#include "correlations.h"
#include "rootoutput.h"

#include "TH1F.h"
#include "TH2S.h"
#include "TFile.h"
#include "TF1.h"

#define XSECPATH "/home/Li6Webb/Desktop/Li6Plus2IAS/li6sim/input/"

using namespace std;

int main(int argc, char *argv[]) {
	
	/**** INPUT ARGUMENTS ****/
  
	string suffix = "";      // output file suffix

	// Total incoming beam energy in MeV, also used for the Fresco simulation.
	// If this changes, make sure to redo the Fresco simulations!
	double Ebeam = 42.82126; // Brho=0.8331 TM (was 42.7 MeV before I revisited Nic's experiment notebooks)

	double Ex     = 2.186; // excitation energy of parent fragment in MeV
	double gamma  = 0.024; // width of excited state of parent fragment in MeV

	// Default physical experiment parameters
	double distanceFromTarget = 235;          // distance of Gobbi from the target in mm
	string suffix             = "_dalpha_3+"; // output file suffix

  // Check for command line arguments, set default values if none are given
	//	- arg. 1 = beam energy in MeV
	//	- arg. 2 = Gobbi distance from target in mm
	//	- arg. 3 = intrinsic state width in MeV
	if (argc >= 3) {
		Ebeam = stod(argv[1]);
		distanceFromTarget = stod(argv[2]);
	}
	else {
		cout << "WARNING: DEFAULT INPUT PARAMETERS BEING USED" << endl;
	}

	// Optional third argument for state width, can supply the first two
	// without this one if desired.
	if (argc == 4)
		gamma = stod(argv[3]);

	// Add simulation parameters to output file suffix, making sure to remove trailing
	// zeros and decimal points.
	string strE = to_string(Ebeam);
	strE.erase(strE.find_last_not_of('0') + 1, string::npos);
	strE.erase(strE.find_last_not_of('.') + 1, string::npos);
	suffix += "_" + strE + "MeV";
	string strDist = to_string(distanceFromTarget);
	strDist.erase(strDist.find_last_not_of('0') + 1, string::npos);
	strDist.erase(strDist.find_last_not_of('.') + 1, string::npos);
	suffix += "_" + strDist + "mm";

	// In case of decimal beam energy, replace '.' with '-' for file prefix purposes
	replace(strE.begin(), strE.end(), '.', '-');

	/**** SETUP AND INITIALIZATION ****/

	bool useRealP = false;  // true means use real angle and energies of fragment
                          // for event reconstruction, to check effect of
                          // detector resolution

	if (useRealP) suffix += "_real";
	if (gamma == 0.) suffix += "_zeroWidth";
	//suffix += "_perfTarg_noResolution2";

	cout << "suffix = " << suffix << " | Ex = " << Ex << " | gamma = " << gamma << endl;

	// Q value here is calculated opposite from how it should to be (according to 
	// Lee's lecture notes), so I added a negative sign to the console output
	double Q = mass_d + mass_alpha - mass_6Li;
  cout << "Q " << -1 * Q << endl;

	// Physical experiment parameters
  float thickness           = 3.026;   // target thickness in mg/cm^2
  float CsiRes              = 0.02;    // resolution of Csi not needed for Si-Si (was 0.00888)
	float b                   = 8.;      // mm beam axis to Gobbi frame dimension,
	float RadiusCollimator    = 0.;      // mm Gobbi collimator outer radius (was 38.1/2.)
  float const targetSize    = 1.0;     // diameter of beam spot size in mm

	// Initialize Gobbi array
	Gobbiarray* gobbi = new Gobbiarray(distanceFromTarget, b, RadiusCollimator);

	// Total cross sections in mb of exit channels for different target excited states from Fresco
	size_t nexits        = 4;                                     // number of exit channels
	double Exts[nexits]  = { 0.0, 3.089443, 3.684507, 3.853807 }; // outgoing target excitation energy for each exit channel

	// Simulation parameters
  int Nevents   = 100000; // events to simulation
	bool einstein = 1;      // switch for newtonian(0) or relativistic(1) kinematics
  float scale   = 1.38;   // scales the magnitude of small angle scattering

	float useRealP_f = (float) useRealP;

	const double Ex_min = 9.5;
  const double dEx = 3.0;

  // Files for energy loss in C target
	string Loss_Li_in_C("Lithium_C.loss");
  string Loss_He_in_C("Helium_C.loss");
  string Loss_d_in_C("Hydrogen_C.loss");

	// Files for energy loss in Si detector
	string Loss_Li_in_Si("Lithium_Si.loss");
	string Loss_He_in_Si("Helium_Si.loss");
	string Loss_d_in_Si("Hydrogen_Si.loss");

	// Initialize fragment objects
	const int Nfrag = 2; // number of decay fragments
  CFrag** frag = new CFrag*[Nfrag];
	frag[0] = new CFrag(1., Mass_d/m0, Loss_d_in_C, Loss_d_in_Si, CsiRes, thickness, gobbi, scale, einstein, useRealP);       // deuteron
	frag[1] = new CFrag(2., Mass_alpha/m0, Loss_He_in_C, Loss_He_in_Si, CsiRes, thickness, gobbi, scale, einstein, useRealP); // alpha

  CFrag *fragBeam = new CFrag(3., Mass_7Li/m0, Loss_Li_in_C, Loss_Li_in_Si, CsiRes, thickness, gobbi, scale, einstein, useRealP);

	// Initialize decay class
  CDecay decay(Nfrag, frag, einstein);

  // Initiallizing the Correlations class reads in the CM cross section from a file
  // and uses that to select a randomized value for phi and theta
	string prefix = "7li12c_e" + strE;
	string Xsecfiles[nexits] = {
		string(XSECPATH) + prefix + "_xsec_2.out",
		string(XSECPATH) + prefix + "_xsec_3.out",
		string(XSECPATH) + prefix + "_xsec_4.out",
		string(XSECPATH) + prefix + "_xsec_5.out"
	};
	string elasXsecfile = string(XSECPATH) + prefix + "_xsec_1.out";
	Correlations* sampler = new Correlations(Xsecfiles, elasXsecfile, Ebeam, Ex, Exts, nexits, Loss_Li_in_C, thickness);

	// Beam momentum and 1.2% MARS acceptance
	double massE = Ebeam + Mass_7Li;
  double pc0 = sqrt((massE*massE) - (Mass_7Li*Mass_7Li));
  //double P_acceptance = .012;
	double P_acceptance = 0.; // use this if not MARS acceptance not relevant

	/**** OUTPUT FILE AND HISTOGRAMS ****/

  RootOutput output(suffix, Nfrag);
	
	/**** MAIN EVENT LOOP ****/

	int Nstuck = 0;
  int Ndet = 0;
  int Nbeamscat = 0;
  for (int index = 0; index < Nevents; index++) {
    // progress updates
    if (index % 1000 == 0) cout << '\xd' <<  index << " of " << Nevents << flush;

		output.Clear();

    // distance in target that produced has to pass to get out
		double rand     = decay.ran.Rndm();
		double inthick  = thickness * rand;
    double outthick = thickness * (1. - rand);

    // beam spot at target
    double rTarget = sqrt(decay.ran.Rndm()) * targetSize / 2.;
    double theta = 2. * pi * decay.ran.Rndm();
    float xTarget = rTarget * cos(theta);
    float yTarget = rTarget * sin(theta);

    // need to re-randomize the angles for each passthrough
    sampler->randomAngles(inthick);
		output.SetSampledValues(&sampler->sampledValues); // save info on beam primary distributions

    //add kinematics to beam
    //simulate MARS having +-1.2% acceptance range
    double pc = pc0 * (1. - P_acceptance / 2.) + decay.ran.Rndm() * P_acceptance * pc0;
    double Ebeam = sqrt((pc * pc) + (Mass_7Li * Mass_7Li)) - Mass_7Li;

		/**** BEAM FRAGMENT PHYSICS ****/

		// set angular properties of beam fragment for elastic scattering case
		double thetaElastic = sampler->sampledValues.GetThetaElasticRad();
		double phi = sampler->sampledValues.GetPhiRad();
    fragBeam->real->SetTheta(thetaElastic);
    fragBeam->real->SetPhi(phi);
    fragBeam->real->SetEnergy(Ebeam);
    fragBeam->real->getVelocity(&einstein); //calculates v, pc & components from energy and angles

    // determine if the beam hits the detector
		fragBeam->targetInteraction(outthick,thickness);
		fragBeam->SiliconInteraction();
		int beamhit = fragBeam->hit(xTarget,yTarget);
		double x, y;
		if (beamhit) {
			x = fragBeam->recon->GetX() / 10.;
			y = fragBeam->recon->GetY() / 10.;
			fragBeam->Egain(thickness * 0.5);
			output.SetElastic(fragBeam->FrontEnergy, fragBeam->DeltaEnergy, fragBeam->recon->GetEnergy(), x, y, fragBeam->recon->GetTheta()*rad_to_deg);
			output.SetIsElasticHit(true);
			Nbeamscat++;
		}

		/**** PARENT FRAGMENT PHYSICS ****/

    // decay parent fragment, add sets velocity vectors of fragments to the seperation
    decay.Mode2Body(Ex, gamma, Q);
		output.SetErelP(decay.ET);

    // transfrom decay vectors to lab frame by adding initial velocity of parent Li6 to all fragments
		double VVparent[3];
    VVparent[0] = sampler->sampledValues.VppX; // x
    VVparent[1] = sampler->sampledValues.VppY; // y
    VVparent[2] = sampler->sampledValues.VppZ; // z
    for (int i = 0; i < Nfrag; i++) frag[i]->AddVelocity(VVparent);

		// Save real Erel post lab frame boost, as sanity check
		output.SetErelPRecon(decay.getErelReal());

		// Save real charged fragment information
		output.SetRealFragment(0, frag[1]->FrontEnergy, frag[1]->DeltaEnergy, frag[1]->real->GetEnergy(), 0., 0., frag[1]->real->GetTheta()*rad_to_deg);
		output.SetRealFragment(1, frag[2]->FrontEnergy, frag[2]->DeltaEnergy, frag[2]->real->GetEnergy(), 0., 0., frag[2]->real->GetTheta()*rad_to_deg);

    // Interaction of fragements in target and silicon detector materials
		// Calculates energy loss in target, change in scatter angle, and
		// wheter fragment is stopped within target
    for (int i = 0; i < Nfrag; i++) {
      frag[i]->targetInteraction(outthick, thickness);
      frag[i]->SiliconInteraction();
    }

		// check for and skip deuterons that punch through back Si layer
		// 20.8911 value is from Lise++ with deuteron and 1.5 mm of Si
    if (frag[0]->FrontEnergy > 20.85) {
			output.Fill();
			continue;
		}

    // detect fragments, and skip if not all fragments detected
    int nhit = 0;
    int ishit = 0;
		int stripx[Nfrag];
    int stripy[Nfrag];
    for (int i = 0; i < Nfrag; i++) {
      ishit = frag[i]->hit(xTarget, yTarget);
			frag[i]->getStripHit(stripx, stripy, i);
      nhit += ishit;
      if (ishit)
        output.DEE->Fill(frag[i]->DeltaEnergy, frag[i]->FrontEnergy);
      if (ishit == -1)
        Nstuck++;
    }

    if (nhit != Nfrag) {
			output.Fill();
			continue;
		}

    // if seperation energy is small, make sure they hit different silicon strips
    // collect what strips are hit
    // loop through all pairs of strips
		bool doublehit = false;
    for (int i = 0; i < Nfrag; i++) {
      for (int j = i + 1; j < Nfrag; j++) {
        // check if it double hit strip
        if (stripx[i] == stripx[j] || stripy[i] == stripy[j])
          doublehit = true; // use to skip after loop
      }

      if (stripx[i] == 32 || stripy[i] == 32) {
        cout << "stripx " << stripx[i] << "   stripy " << stripy[i] << endl;
        doublehit = true;
      }
    }

    if (doublehit) {
			output.Fill();
			continue;
		}

    // We have a detection
		output.SetIsFragDet(true);
    Ndet++;

    for (int i = 0; i < Nfrag; i++)
      frag[i]->Egain(thickness * 0.5);

		output.alphaenergy->Fill(frag[1]->recon->GetEnergy());
    output.protonenergy->Fill(frag[0]->recon->GetEnergy());

    // Get reconstructed  relative energy between fragements
    float Erel_S = useRealP_f ? decay.getErelReal() : decay.getErelRecon();

    //get reconstructed excitation energy
    float Ex_S = Erel_S + Q;

		decay.plfRecon->SetEnergy(Erel_S);
		decay.plfRecon->RadToDeg();
		output.SetEx(Ex_S);
		output.SetCosThetaH(decay.cos_thetaH);
    output.hist_Erel_thetaH->Fill(Erel_S, decay.cos_thetaH);

    //look at transverse emisson for better resolutions
    if (fabs(decay.cos_thetaH) < 0.7) output.hist_Ex_trans->Fill(Ex_S);
    if (fabs(decay.cos_thetaH) < 0.5) output.hist_Ex_trans_narrow->Fill(Ex_S);

    output.hist_Ex_DE->Fill(Ex_S, frag[1]->FrontEnergy);
		output.SetReconValues(decay.plfRecon->GetKinematicValues());

		x = frag[0]->recon->GetX()/10.;
    y = frag[0]->recon->GetY()/10.;
    output.protonXY_S->Fill(x,y);
		output.SetReconFragment(0, frag[0]->FrontEnergy, frag[0]->DeltaEnergy, frag[0]->recon->GetEnergy(), x, y, frag[0]->recon->GetTheta()*rad_to_deg);

    x = frag[1]->recon->GetX()/10.;
    y = frag[1]->recon->GetY()/10.;
    output.coreXY_S->Fill(x,y);
		output.SetReconFragment(1, frag[1]->FrontEnergy, frag[1]->DeltaEnergy, frag[1]->recon->GetEnergy(), x, y, frag[1]->recon->GetTheta()*rad_to_deg);

		output.Fill();
  }

  cout << endl;
  cout << endl;

  cout << "(d + alpha) coincidence efficiency = " << (float)Ndet/(float)Nevents << endl;
  cout << "(d + alpha) was stuck in target = " << Nstuck << ", fraction = " << (float)Nstuck/(float)Nevents << endl;
  cout << "7Li beam elastic scatter det = " << (float)Nbeamscat/(float)Nevents << endl;
  
  //make fit to measure resolution of Invarient mass peak
  //TF1 * fit = new TF1("fit", "gaus", Ex_min, Ex_min+dEx);
  //hist_Ex_S->Fit("fit", "R");
  //double mean = fit->GetParameter(1);
  //double sigma = fit->GetParameter(2);
  //cout << "mean " << mean << " sigma " << sigma << endl;

	//clean up, clean up
  delete[] frag;
  //everybody everywhere
  delete sampler;
  //clean up, clean up
	delete fragBeam;
	//everybody do your share
	delete gobbi;

  //beep at me when finished (sadly doesn't work anymore)
  //cout << "\a";
}
