#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <sstream>
#include <string>
#include <iomanip>
#include "cuts.h"
#include "Gobbiarray.h"
#include "frag.h"
#include "decay.h"
#include "constants.h"
#include "correlations.h"
#include "rootoutput.h"

#include <Math/MinimizerOptions.h>
#include <TH1F.h>
#include <TH2S.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TInterpreter.h>

using namespace std;

int main(int argc, char *argv[]) {

	// Make sure dictionary is properly linked and stuff
	cout << "Loading simlib shared library from: " << string(SOFILE) << endl;
	cout << "Loading simlib ROOT dictionary from: " << string(PCMFILE) << endl;
	gSystem->Load(SOFILE);
	TInterpreter::Instance()->AddIncludePath(PCMFILE);

	// Default minimizer
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Migrad");
	
	/**** INPUT ARGUMENTS ****/

	// Total incoming beam energy in MeV, also used for the Fresco simulation.
	// If this changes, make sure to redo the Fresco simulations!
	//double Ebeam = 42.82126; // Brho=0.8331 TM (was 42.7 MeV before I revisited Nic's experiment notebooks)
	double Ebeam = 56;

	double Ex    = 2.186; // excitation energy of parent fragment in MeV
	double gamma = 0.024; // width of excited state of parent fragment in MeV

	// Default physical experiment parameters
	double distanceFromTarget = 90;           // distance of Gobbi from the target in mm, 235 for Nic's experiment
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

	// Q value here is calculated opposite from how it should to be (according to 
	// Lee's lecture notes), so I added a negative sign to the console output
	double Q = mass_d + mass_alpha - mass_6Li;
	cout << "Q " << -1 * Q << endl;

	// Physical experiment parameters
	float thickness           = 17.575;  // target thickness in mg/cm^2 (3.026 for Nic's experiment)
	float density             = 3515.;   // target density in mg/cm^3 (2260 for graphite, 3515 for diamond)
	float CsiRes              = 0.02;    // resolution of Csi not needed for Si-Si (was 0.00888)
	float b                   = 8.;      // mm beam axis to Gobbi frame dimension,
	float RadiusCollimator    = 0.;      // mm Gobbi collimator outer radius (was 38.1/2.)
	float const targetSize    = 1.0;     // diameter of beam spot size in mm
	float quenchRecoil = .0;             // quenching factor of recoil in target (min 0 = full energy deposition, max 1 = no energy deposition)
	
	float diamondResFWHM = 0.1;                                     // diamond detector energy resolution (FWHM) (MeV)
	float diamondRes = diamondResFWHM * .5f / sqrt(2.f * log(2.f)); // diamond detector energy resolution (sigma) (MeV)

	// Initialize Gobbi array
	shared_ptr<Gobbiarray> gobbi = make_shared<Gobbiarray>(distanceFromTarget, thickness / density * 10., b, RadiusCollimator);

	// Total cross sections in mb of exit channels for different target excited states from Fresco
	size_t nexits       = 4;                                     // number of exit channels
	vector<double> Exts = { 0.0, 3.089443, 3.684507, 3.853807 }; // outgoing target excitation energy for each exit channel

	// Simulation parameters
	int Nevents   = 100000; // events to simulation
	bool einstein = 1;      // switch for newtonian(0) or relativistic(1) kinematics
	float scale   = 1.38;   // scales the magnitude of small angle scattering
	
	// List of target reaction z position test points (in % from upstream side, must range 0-1)
	size_t nTestPoints{5};
	std::vector<double> thickTests; // MUST BE IN mg/cm^2 (this is automatically filled in the Init function using the set number of test points)
	std::vector<double> dETests;
	
	// Set up target total energy loss function test point vectors
	thickTests.resize(nTestPoints);
	dETests.resize(nTestPoints);
	double h = thickness / (double)(nTestPoints - 1);
	for (int i = 0; i < nTestPoints; i++)
		thickTests[i] = h * (double)i;

	float useRealP_f = (float) useRealP;

	//const double Ex_min = 9.5;
	//const double dEx = 3.0;
	
	/******** FILE SUFFIX CREATION ********/

	ostringstream oss;

	// Source - https://stackoverflow.com/a/46424921
	// Posted by Ivan Folgueira Bande, modified by community. See post 'Timeline' for change history
	// Retrieved 2026-03-30, License - CC BY-SA 4.0
	oss << noshowpoint << (diamondResFWHM * 1000.);
	suffix += "_diamond" + oss.str() + "keV";
	
	suffix += "_loss0-005step";
	
	// Source - https://stackoverflow.com/a/7623670
	// Posted by Darcy Rayner
	// Retrieved 2026-03-30, License - CC BY-SA 3.0
	oss.str("");
	oss.clear(); // Clear state flags.

	oss << noshowpoint << quenchRecoil;
	suffix += "_" + oss.str() + "recoilquench";

	// Source - https://stackoverflow.com/a/2896627
	// Posted by Kirill V. Lyadvinsky, modified by community. See post 'Timeline' for change history
	// Retrieved 2026-03-30, License - CC BY-SA 3.0
	replace(suffix.begin(), suffix.end(), '.', '-');
	
	if (useRealP) suffix += "_real";
	if (gamma == 0.) suffix += "_zeroWidth";
	//suffix += "_perfTarg_noResolution2";

	cout << "suffix = " << suffix << " | Ex = " << Ex << " | gamma = " << gamma << endl;

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
	vector<shared_ptr<CFrag>> frag;
	frag.resize(Nfrag);
	frag[0] = make_shared<CFrag>(1., Mass_d/m0, Loss_d_in_C, Loss_d_in_Si, CsiRes, thickness, gobbi, scale, einstein, useRealP);       // deuteron
	frag[1] = make_shared<CFrag>(2., Mass_alpha/m0, Loss_He_in_C, Loss_He_in_Si, CsiRes, thickness, gobbi, scale, einstein, useRealP); // alpha

	CFrag *fragBeam = new CFrag(3., Mass_7Li/m0, Loss_Li_in_C, Loss_Li_in_Si, CsiRes, thickness, gobbi, scale, einstein, useRealP);

	// Initialize decay class
	CDecay decay(Nfrag, frag, einstein);

	// Initiallizing the Correlations class reads in the CM cross section from a file
	// and uses that to select a randomized value for phi and theta
	string prefix = "7li12c_e" + strE;
	vector<string> Xsecfiles = {
		string(XSECPATH) + "/" + prefix + "_xsec_2.out",
		string(XSECPATH) + "/" + prefix + "_xsec_3.out",
		string(XSECPATH) + "/" + prefix + "_xsec_4.out",
		string(XSECPATH) + "/" + prefix + "_xsec_5.out"
	};
	string elasXsecfile = string(XSECPATH) + "/" + prefix + "_xsec_1.out";
	Correlations* sampler = new Correlations(Xsecfiles, elasXsecfile, Ebeam, Ex, Exts, nexits, Loss_Li_in_C, thickness);

	// Beam momentum and 1.2% MARS acceptance
	double massE = Ebeam + Mass_7Li;
	double pc0 = sqrt((massE*massE) - (Mass_7Li*Mass_7Li));
	//double P_acceptance = .012;
	double P_acceptance = 0.; // use this if not MARS acceptance not relevant
	
	// Helper function to calculate the total target energy loss using a given target reaction z position, the
	// reconstructed Gobbi variables, and other known experimental parameters. Here, the only unknown value
	// is the target reaction z position. The idea is to calculate this for a few different points, fit a
	// function, and then use this fit function to invert the non-trivial target energy loss function. THIS IS 
	// THE FULLY QUENCHED VERSION.
	//     inthick - reaction z position in target from upstream side (mg / cm^2).
	auto CalcTargELossQuenched = [&](double inthick) {
		double dETarg   = Ebeam - fragBeam->loss_C->getEout(Ebeam, inthick);
		double outthick = thickness - inthick;
		CFrame* alphFrame = frag[1]->recon;
		CFrame* deutFrame = frag[0]->recon;
		if (useRealP) {
			alphFrame = frag[1]->real;
			deutFrame = frag[0]->real;
		}
		double dEAlpha  = frag[1]->loss_C->getEin(alphFrame->GetEnergy(), outthick / cos(alphFrame->GetTheta())) - alphFrame->GetEnergy();
		double dEDeut = frag[0]->loss_C->getEin(deutFrame->GetEnergy(), outthick / cos(deutFrame->GetTheta())) - deutFrame->GetEnergy();
		return dETarg + dEAlpha + dEDeut;
	};

	// Helper function to calculate the total target energy loss using a given target reaction z position, the
	// reconstructed Gobbi variables, and other known experimental parameters. Here, the two unknown values
	// are the reaction depth and the recoil excitation energy. The latter is determined externally per event
	// by gating on diamond energy and total Gobbi energy. This function is then applied to a few different 
	// test reaction depths. The results are fit with a function, which is then inverted to convert a diamond
	// energy into a reaction depth. THIS VERSION HAS ZERO QUENCHING.
	//     inthick - reaction z position in target from upstream side (mg / cm^2)
	//     ext     - excitation energy of target recoil nucleus, determined externally
	auto CalcTargELoss = [&](double inthick, double ext) {
		double outthick = thickness - inthick;
		CFrame* temp_frags[Nfrag];
		double temp_dEtarg = Ebeam;
		for (size_t i = 0; i < Nfrag; i++) {
			temp_frags[i] = new CFrame(*(useRealP ? frag[i]->real : frag[i]->recon));
			temp_dEtarg -= temp_frags[i]->GetEnergy();
			temp_dEtarg += frag[i]->EgainHelper(outthick / cos(temp_frags[i]->GetTheta()), temp_frags[i]);
		}

		temp_dEtarg -= decay.getErel(temp_frags) + sampler->GetQValue() - ext;

		return temp_dEtarg;
	};

	/**** OUTPUT FILE AND HISTOGRAMS ****/

	RootOutput output(suffix, Nfrag, false);
	
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

		//add kinematics to beam
		//simulate MARS having +-1.2% acceptance range
		if (P_acceptance != 0.) {
			double pc = pc0 * (1. - P_acceptance / 2.) + decay.ran.Rndm() * P_acceptance * pc0;
			Ebeam = sqrt((pc * pc) + (Mass_7Li * Mass_7Li)) - Mass_7Li;
		}

		// Calculate energy dropped by beam in target, save for later
		double beamEloss = Ebeam - fragBeam->loss_C->getEout(Ebeam, inthick);
		double dETarg = beamEloss;

		// beam spot at target
		double rTarget = sqrt(decay.ran.Rndm()) * targetSize / 2.;
		double theta = 2. * pi * decay.ran.Rndm();
		float xTarget = rTarget * cos(theta);
		float yTarget = rTarget * sin(theta);

		// need to re-randomize the angles for each passthrough
		sampler->randomAngles(inthick);
		output.SetSampledValues(&sampler->sampledValues); // save info on beam primary distributions

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
		int beamhit = fragBeam->hit(inthick, xTarget, yTarget);
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
		output.SetRealFragment(0, frag[0]->FrontEnergy, frag[0]->DeltaEnergy, frag[0]->real->GetEnergy(), 0., 0., frag[0]->real->GetTheta()*rad_to_deg);
		output.SetRealFragment(1, frag[1]->FrontEnergy, frag[1]->DeltaEnergy, frag[1]->real->GetEnergy(), 0., 0., frag[1]->real->GetTheta()*rad_to_deg);

		// Interaction of fragements in target and silicon detector materials
		// Calculates energy loss in target, change in scatter angle, and
		// wheter fragment is stopped within target
		for (int i = 0; i < Nfrag; i++) {
			dETarg += frag[i]->real->GetEnergy();
			frag[i]->targetInteraction(outthick, thickness);
			dETarg -= frag[i]->real->GetEnergy();
			frag[i]->SiliconInteraction();
		}

		// Add recoil energy loss to total energy loss in target, assuming zero quenching and all energy is deposited
		// Here, the recoil quenching factor is how much quenching there is. In other words, a factor of 0 means
		// that there is no quenching, so the full recoil energy is depositied as measurable light. A factor of 1
		// means that there is 100% quenching, so the recoil deposits no measurable light. Hence, why the recoil
		// factor is subtracted from 1 here.
		dETarg += (1. - quenchRecoil) * sampler->sampledValues.recoilEnergy;

		// Reconstruct reaction position in target from total target energy loss
		double dETargRecon = dETarg + decay.ran.Gaus(0., diamondRes);
		//double inthickreconavg = min(max((dETargRecon - 5.10193) / 0.184146, 0.), (double)thickness); // linear function from fitting dETarg vs. inthick
		double inthickreconavg = -1; // average linear formula is probably not useful or valid here

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
			ishit = frag[i]->hit(inthick, xTarget, yTarget);
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

		/******** REACTION DEPTH RECONSTRUCTION ********/

		// First, identify the target recoil nucleus' excitation energy
		double sumGobbiE = 0.;
		for (size_t i = 0; i < Nfrag; i++)
			sumGobbiE += (useRealP ? frag[i]->real : frag[i]->recon)->GetEnergy();
		double Ext = -1;
		if (cut_gs.IsInside(sumGobbiE, dETargRecon)) {
			Ext = Exts[0];
			output.SetValidExt(true);
		}
		else if (cut_52p.IsInside(sumGobbiE, dETargRecon)) {
			Ext = Exts[3];
			output.SetValidExt(true);
		}

		// Calculate target total energy loss test points before addback
		dETests.clear();
		dETests.resize(thickTests.size());

		for (int i = 0; i < thickTests.size(); i++)
			dETests[i] = CalcTargELoss(thickTests[i], Ext);
		output.SetDETests(dETests);

		// Quadratic fit to target total energy loss function
		TF1 elossFit("quadratic", "([0]*x*x)+([1]*x)+[2]");
		elossFit.SetRange(-0.1, thickness+0.1);
		elossFit.SetParameter(0, 1);
		elossFit.SetParLimits(0, 0, 1);
		
		TGraph elossGraph(thickTests.size(), &thickTests[0], &dETests[0]);
		elossGraph.Fit(&elossFit, "NRQ", "", -0.1, thickness+0.1);

		// Invert target energy loss function and calculate target reaction position
		double a = elossFit.GetParameter(0);
		double b = elossFit.GetParameter(1);
		double c = elossFit.GetParameter(2) - dETargRecon;
		double disc = (b*b) - (4. * a * c);
		double inthickrecon = -1;
		if (disc >= 0) {
			inthickrecon = 0.5 * (-b + sqrt(disc)) / a;
			inthickrecon = min(max(inthickrecon, 0.), (double)thickness); // clamp value to within target dimensions
		}
		else cout << "Discriminant < 0, something went very wrong!" << endl;
		//inthickrecon = inthick; // uncomment this if you want a perfect reaction position reconstruction
		output.SetTargetEloss(dETarg, dETargRecon, inthick, inthickreconavg, inthickrecon);

		// Energy addback for target
		double fragEgain = 0.;
		for (int i = 0; i < Nfrag; i++) {
			fragEgain -= (useRealP ? frag[i]->real : frag[i]->recon)->GetEnergy();
			fragEgain += frag[i]->Egain((thickness - ((disc >= 0) ? inthickrecon : inthickreconavg)) / cos((useRealP ? frag[i]->real : frag[i]->recon)->GetTheta()));
			//frag[i]->Egain(thickness * 0.5);
		}

		// Calculate Eloss-Egain for the fragments in the target, for debugging purposes
		dETarg -= beamEloss;
		output.SetFragElossEgain(dETarg, fragEgain);

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
	//everybody everywhere
	delete sampler;
	//clean up, clean up
	delete fragBeam;
	//everybody do your share

	//beep at me when finished (sadly doesn't work anymore)
	//cout << "\a";
}
