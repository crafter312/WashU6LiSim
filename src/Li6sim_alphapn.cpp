#include "Li6sim_alphapn.h"

#include "constants.h"

#include <iostream>

using namespace std;

Li6sim_alphapn::Li6sim_alphapn(double Eb, double GobbiDist, double _ex, double width, string suf)
	: Ebeam(Eb), distanceFromTarget(GobbiDist), Ex(_ex), gamma(width), suffix(suf) {

	gobbi = make_shared<Gobbiarray>(GobbiDist, b, RadiusCollimator);

	if(suffix != "") suffix = "_" + suffix;

	// Add simulation parameters to output file suffix, making sure to remove
	// trailing zeros and decimal points.
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

	if (gamma == 0.) suffix += "_zeroWidth";

	// Q value here is calculated opposite from how it should to be (according to 
	// Lee's lecture notes), so I added a negative sign to the console output
	Q = mass_alpha + mass_p + mass_n - mass_6Li;

	// Set exit channel info
	nexits = 4;
	Exts = { 0.0, 3.089443, 3.684507, 3.853807 };

	// Initiallizing the Correlations class reads in the CM cross section from a file
	// and uses that to select a randomized value for phi and theta
	string prefix = "7li12c_e" + strE;
	Xsecfiles = {
		string(XSECPATH) + "/" + prefix + "_xsec_2.out",
		string(XSECPATH) + "/" + prefix + "_xsec_3.out",
		string(XSECPATH) + "/" + prefix + "_xsec_4.out",
		string(XSECPATH) + "/" + prefix + "_xsec_5.out"
	};
	elasXsecfile = string(XSECPATH) + "/" + prefix + "_xsec_1.out";

	// Set number of fragments
	Nfrag = 3;

	// Beam momentum
	double massE = Ebeam + Mass_7Li;
	pc0 = sqrt((massE*massE) - (Mass_7Li*Mass_7Li));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Li6sim_alphapn::Li6sim_alphapn(double Eb, double GobbiDist, double _ex, double width, string suf, float b0, float radColl) {
	b = b0;
	RadiusCollimator = radColl;
	Li6sim_alphapn(Eb, GobbiDist, _ex, width, suf);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Li6sim_alphapn::~Li6sim_alphapn() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

string Li6sim_alphapn::Init() {

	// First, double check the exit channel vector lengths
	if((Exts.size() != nexits) || (Xsecfiles.size() != nexits))
		return "Exit channel vector lengths don't match that expected of nexits = " + to_string(nexits);

	if(useRealP) suffix += "_real";

	// It's important here that the heavy fragment (alpha) is the last fragment, and
	// also that the neutron is the first fragment so that I can skip it.
	frag = {
		make_shared<CFrag>(0., Mass_n/m0, Loss_n_in_C, Loss_n_in_Si, GobbiRes, thickness, gobbi, scale, einstein, useRealP),      // neutron
		make_shared<CFrag>(1., Mass_p/m0, Loss_p_in_C, Loss_p_in_Si, GobbiRes, thickness, gobbi, scale, einstein, useRealP),      // proton
		make_shared<CFrag>(2., Mass_alpha/m0, Loss_He_in_C, Loss_He_in_Si, GobbiRes, thickness, gobbi, scale, einstein, useRealP) // alpha
	};

	// Initialize beam fragment
	fragBeam = make_unique<CFrag>(3., Mass_7Li/m0, Loss_Li_in_C, Loss_Li_in_Si, GobbiRes, thickness, gobbi, scale, einstein, useRealP);

	// Initialize decay class
	decay = make_unique<CDecay>(Nfrag, frag, einstein);

	// Initialize sampler class
	sampler = make_unique<Correlations>(Xsecfiles, elasXsecfile, Ebeam, Ex, Exts, nexits, Loss_Li_in_C, thickness);

	// Update momentum acceptance, in case 1.2% MARS acceptance was enabled
	P_acceptance = mars * .012;

	return "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Li6sim_alphapn::PrintSettings() {
	cout << "=================================================================" << endl;
	cout << "Charged particle simulation settings:" << endl;
	cout << "suffix = " << suffix << " | Ex = " << Ex << " | gamma = " << gamma << endl;

	// The Q value used in the simulation is calculated opposite from how it
	// should to be (according to Lee's lecture notes), so I added a negative
	// sign to the console output.
	cout << "Q = " << -1 * Q << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

string Li6sim_alphapn::DoSingleEventPreNeutron(RootOutput& output) {

	// Make sure important variables are correctly initialized
	if(frag.size() != Nfrag)
		return "Fragment vector size doesn't match that expected of Nfrag = " + to_string(Nfrag);
	for(int i = 0; i < Nfrag; i++)
		if(!frag[i]) return "Fragment pointer at index " + to_string(i) + " not initialized properly";
	if(!fragBeam)
		return "Beam fragment pointer not initialized properly";
	if(!decay)
		return "Decay class pointer not initialized properly";
	if(!sampler)
		return "Sampler class pointer not initialized properly";

	output.Clear();

	// Distance in target that produced has to pass to get out
	double rand		  = decay->ran.Rndm();
	double inthick	= thickness * rand;
	double outthick = thickness * (1. - rand);

	// Beam spot at target
	double rTarget = sqrt(decay->ran.Rndm())*targetSize/2.;
	double theta = 2.*acos(-1.)*decay->ran.Rndm();
	float xTarget = rTarget*cos(theta);
	float yTarget = rTarget*sin(theta);

	// Need to re-randomize the angles for each passthrough
	sampler->randomAngles(inthick);
	output.SetSampledValues(&sampler->sampledValues); // save info on beam primary distributions

	// Simulate MARS having +-1.2% acceptance range, if enabled
	double pc = pc0*(1.-P_acceptance/2.) + decay->ran.Rndm()*P_acceptance*pc0;
	double Ebeam = sqrt(pow(pc,2) + pow(Mass_7Li,2)) - Mass_7Li;

	/**** BEAM FRAGMENT PHYSICS ****/

	// Set angular properties of beam fragment for elastic scattering case
	double thetaElastic = sampler->sampledValues.GetThetaElasticRad();
	double phi = sampler->sampledValues.GetPhiRad();
	fragBeam->real->SetTheta(thetaElastic);
	fragBeam->real->SetPhi(phi);
	fragBeam->real->SetEnergy(Ebeam);
	fragBeam->real->getVelocity(&einstein); //calculates v, pc & components from energy and angles

	// Determine if the beam hits the detector
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

	// Decay parent fragment, add sets velocity vectors of fragments to the seperation
	decay->ModeMicroCanonical(Ex, gamma, Q);
	output.SetErelP(decay->ET);

	// Transfrom decay vectors to lab frame by adding initial velocity of parent Li7 to all fragments
	double VVparent[3];
	VVparent[0] = sampler->sampledValues.VppX; // x
	VVparent[1] = sampler->sampledValues.VppY; // y
	VVparent[2] = sampler->sampledValues.VppZ; // z
	for (int i = 0; i < Nfrag; i++) frag[i]->AddVelocity(VVparent);

	// Save real Erel post lab frame boost, as sanity check
	output.SetErelPRecon(decay->getErelReal());

	// Save real charged fragment information
	output.SetRealFragment(0, frag[1]->FrontEnergy, frag[1]->DeltaEnergy, frag[1]->real->GetEnergy(), 0., 0., frag[1]->real->GetTheta()*rad_to_deg);
	output.SetRealFragment(1, frag[2]->FrontEnergy, frag[2]->DeltaEnergy, frag[2]->real->GetEnergy(), 0., 0., frag[2]->real->GetTheta()*rad_to_deg);

	/**** CHARGED FRAGMENT RECONSTRUCTION ****/

	// Interaction of fragements in target and silicon detector materials
	// Calculates energy loss in target, change in scatter angle, and
	// wheter fragment is stopped within target
	for (int i = 1; i < Nfrag; i++) {
		frag[i]->targetInteraction(outthick, thickness);
		frag[i]->SiliconInteraction();
	}

	// check for and skip protons that punch through back Si layer
	// 15.5 value is from Lise++ with proton and 1.5 mm of Si
	if (frag[1]->FrontEnergy > 15.5) {
		output.Fill();
		Npunch++;
		return "";
	}

	// detect fragments, and skip if not all fragments detected
	int nhit = 0;
	int ishit = 0;
	int stripx[Nfrag - 1];
	int stripy[Nfrag - 1];
	for (int i = 1; i < Nfrag; i++) {
		ishit = frag[i]->hit(xTarget, yTarget);
		frag[i]->getStripHit(stripx, stripy, i);
		nhit += ishit;
		if (ishit)
			output.DEE->Fill(frag[i]->DeltaEnergy, frag[i]->FrontEnergy);
		if (ishit == -1)
			Nstuck++;
	}

	if (nhit != Nfrag - 1) {
		output.Fill();
		Nmiss++;
		return "";
	}

	// if seperation energy is small, make sure they hit different silicon strips
	// collect what strips are hit
	// loop through all pairs of strips
	bool doublehit = false;
	for (int i = 1; i < Nfrag; i++) {
		for (int j = i + 1; j < Nfrag; j++) {
			// check if it double hit strip
			if (stripx[i] == stripx[j] || stripy[i] == stripy[j])
				doublehit = true; // use to skip after loop
		}

		if (stripx[i] == 32 || stripy[i] == 32) {
			cout << "stripx " << stripx[i] << "	 stripy " << stripy[i] << endl;
			doublehit = true;
		}
	}

	if (doublehit) {
		output.Fill();
		N2Hit++;
		return "";
	}
	
	// We have a detection
	output.SetIsFragDet(true);
	Ndet++;

	// Energy addback for half target
	for (int i = 1; i < Nfrag; i++) {
		frag[i]->Egain(thickness * 0.5);
	}

	// Output of charged fragment information
	output.protonenergy->Fill(frag[1]->recon->GetEnergy());
	output.alphaenergy->Fill(frag[2]->recon->GetEnergy());

	x = frag[1]->recon->GetX() / 10.;
	y = frag[1]->recon->GetY() / 10.;
	output.protonXY_S->Fill(x,y);
	output.SetReconFragment(0, frag[1]->FrontEnergy, frag[1]->DeltaEnergy, frag[1]->recon->GetEnergy(), x, y, frag[1]->recon->GetTheta()*rad_to_deg);

	x = frag[2]->recon->GetX() / 10.;
	y = frag[2]->recon->GetY() / 10.;
	output.coreXY_S->Fill(x,y);
	output.SetReconFragment(1, frag[2]->FrontEnergy, frag[2]->DeltaEnergy, frag[2]->recon->GetEnergy(), x, y, frag[2]->recon->GetTheta()*rad_to_deg);

	return "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

string Li6sim_alphapn::DoSingleEventPostNeutron(RootOutput& output) {

	/**** NEUTRON RECONSTRUCTION ****/

	// Detector geometry not implemented yet, assume distance of 1 m
	double neutDist = 100; // cm

	// Fold in 1ns timing resolution
	double neutV = frag[0]->real->GetVelocity(); // cm/ns
	double neutT = neutDist / neutV; // ns
	if (!useRealP) neutT += decay->ran.Gaus(0., neutTRes);
	output.SetTNeut(neutT);
	frag[0]->recon->SetVelocity(neutDist / neutT);

	// Assume other values are exact
	frag[0]->recon->SetTheta(frag[0]->real->GetTheta());
	frag[0]->recon->SetPhi(frag[0]->real->GetPhi());
	frag[0]->recon->Sph2CartV();
	frag[0]->recon->getEnergy(&einstein);

	//get reconstructed relative energy between fragements
	float Erel_S = useRealP ? decay->getErelReal() : decay->getErelRecon();

	//get reconstructed excitation energy
	float Ex_S = Erel_S + Q;

	decay->plfRecon->SetEnergy(Erel_S);
	decay->plfRecon->RadToDeg();
	output.SetEx(Ex_S);
	output.SetCosThetaH(decay->cos_thetaH);
	output.hist_Erel_thetaH->Fill(Erel_S, decay->cos_thetaH);

	//look at transverse emisson for better resolutions
	if (fabs(decay->cos_thetaH) < 0.7) output.hist_Ex_trans->Fill(Ex_S);
	if (fabs(decay->cos_thetaH) < 0.5) output.hist_Ex_trans_narrow->Fill(Ex_S);

	output.hist_Ex_DE->Fill(Ex_S, frag[2]->FrontEnergy);
	output.SetReconValues(decay->plfRecon->GetKinematicValues());

	output.SetENeut(frag[0]->recon->GetEnergy());
	output.SetThetaNeut(frag[0]->recon->GetTheta()*rad_to_deg);
	output.Fill();

	return "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Li6sim_alphapn::DoFinalThings(int Nevents) {
	cout << endl;
	cout << endl;

	cout << "Skipped event counts:" << endl;
	cout << "p punch throughs = " << (float)Npunch/(float)Nevents << endl;
	cout << "geometry cuts = " << (float)Nmiss/(float)Nevents << endl;
	cout << "double hits = " << (float)N2Hit/(float)Nevents << endl;

	cout << endl;

	cout << "(alpha + p + n) coincidence efficiency = " << (float)Ndet/(float)Nevents << endl;
	cout << "(alpha + p + n) was stuck in target = " << Nstuck << ", fraction = " << (float)Nstuck/(float)Nevents << endl;
	cout << "7Li beam elastic scatter det = " << (float)Nbeamscat/(float)Nevents << endl;
}



