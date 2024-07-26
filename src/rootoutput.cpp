/**
 * Author: Henry Webb (h.s.webb@wustl.edu)
 * Created: 22 May 2024
 */

#include "rootoutput.h"

// Input:
//	nFrags -- number of decay fragments
RootOutput::RootOutput(string suffix, int n) {
	nFrags = n;

	// Define output file
  string fileout = "sim" + suffix + ".root";
  file = new TFile(fileout.c_str(), "RECREATE");
	file->cd();

	/**** Initialize tree ****/

	t = new TTree("t", "t");
  t->Branch("parentSecondary", &parentSecondary);
  t->Branch("thetaNeut", &thetaNeut);
	t->Branch("ErelP", &ErelP);
	t->Branch("ErelS", &ErelS);
	t->Branch("Ex", &Ex);
	t->Branch("cosThetaH", &cosThetaH);
	t->Branch("isElasticHit", &isElasticHit);
	t->Branch("isFragDet", &isFragDet);

	t->Branch("sampler", &sampler);

	// Decay fragments
	chargedFragments = new CFragS[nFrags];
	for (int i=0; i<nFrags; i++)
		chargedFragments[i] = CFragS();
	t->Branch("chargedFragments", &chargedFragments);

	/**** Initialize histograms ****/

	// Primary distributions for parent fragment
  hist_vel_P       = new TH1F("vel_P","vel",100,1.5,4);
  hist_theta_P     = new TH1F("theta_P","theta",200,0,30);
  hist_phi_P       = new TH1F("phi_P","phi",180,0,360);
  kinematic_circle = new TH2F("kinematic_circle", "kin_circle", 200, -1, 1, 120, 2, 3.8);

	// Primary elastic scattering distribution
  hist_theta_beam_P = new TH1F("theta_beam_P","theta",200,0,30);

  // Reconstructed seconday distributions of parent;
  hist_vel_S   = new TH1F("vel_S","vel",100,1.5,4);
  hist_theta_S = new TH1F("theta_S","theta",200,0,30);
  hist_phi_S   = new TH1F("phi_S","phi",180,0,360);

	// Secondary real and reconstructed elastic scattering distributions
  hist_theta_beam_S_sharp = new TH1F("hist_theta_beam_S_sharp","theta",200,0,30);
  hist_theta_beam_S_recon = new TH1F("hist_theta_beam_S_recon","theta",200,0,30);

	// Reconstructed fragment energy distributions
  DEE          = new TH2I("DEE","",800,0,22,500,0,80); // E is x, DE is y
  protonenergy = new TH1F("protonenergy","", 100,0,30);
  alphaenergy  = new TH1F("alphaenergy","", 100,0,40);

	// Angular distribution of heavy fragment
	// -1 or 1 is longitudinal decay, 0 is transverse decay
  cos_thetaH       = new TH1F("cos_thetaH","",100,-1,1);
  hist_Erel_thetaH = new TH2F("Erel_thetaH","",200,0,8,25,-1,1);

	// Decay and excitation energy
  hist_Erel_P          = new TH1I("hist_Erel_P","",200,0,8);
  hist_Erel            = new TH1I("Erel","",200,0,8);
  hist_Ex              = new TH1I("Ex","",400,0,18);
  hist_Ex_trans        = new TH1I("Ex_trans","",400,0,18);
  hist_Ex_trans_narrow = new TH1I("Ex_trans_narrow","",400,0,18);
	hist_Ex_DE           = new TH2I("Ex_DE","", 400,0,18, 100,0,16);

  // Maps of x-y positions of detected fragments
  protonXY_S = new TH2S("protonXY_S","protonxy",100,-10,10,100,-10,10);
  coreXY_S   = new TH2S("coreXY_S","alphaxy",100,-10,10,100,-10,10);

	// Angle of neutron fragment
	hist_neut_theta = new TH1F("neut_theta","",200,0,90);
}

RootOutput::~RootOutput() {
	file->Write();
	file->Close();

	delete[] chargedFragments;
}

void RootOutput::Fill() {
	t->Fill();
}

void RootOutput::Clear() {
	parentSecondary.clear();
	for (int i = 0; i < nFrags; i++)
		chargedFragments[i].clear();
	thetaNeut = -1;
	thetaElastS = -1;
	ErelP = -1;
	ErelS = -1;
	Ex = -1;
	cosThetaH = NAN;
	isElasticHit = -1;
	isFragDet = false;

	sampler.Clear();
}

// Input:
//	v -- total velocity
//	vx -- x velocity
//	vy -- y velocity
//	vz -- z velocity
// 	p -- phi (azimuthal angle) in degrees
//	th -- theta (polar angle) in degrees
// The first six inputs are for the case of inelastic scattering,
// with this instance being for the "primary" distribution (prior to
// detection, reconstruction, etc.)

// Input: SEE PREVIOUS FUNCTION
// These parameters are the same as the previous function, except that
// they are for the "secondary" distribution (after detection,
// reconstruction, etc.).
// Also note that both "p" and "th" should be in radians here!
void RootOutput::SetSecondary(double v, double p, double th) {
	parentSecondary.vel   = v;
	parentSecondary.velx  = v*sin(th)*cos(p);
	parentSecondary.vely  = v*sin(th)*sin(p);
	parentSecondary.velz  = v*cos(th);
	parentSecondary.phi   = p*rad_to_deg;
	parentSecondary.theta = th*rad_to_deg;

	hist_vel_S->Fill(parentSecondary.vel);
  hist_theta_S->Fill(parentSecondary.theta);
  hist_phi_S->Fill(parentSecondary.phi);
}

// Input:
//	n -- the number of charged fragments (not including neutrons)
//	de -- the energy lost in the thin front Si detector
//	e -- the remaining energy lost in the thick back Si detector
//	recE -- reconstructed fragment energy after detection
//	_x -- x position of fragment in detector
//	_y -- y position of fragment in detector
void RootOutput::SetFragment(int n, double de, double e, double recE, double _x, double _y) {
	chargedFragments[n].DE          = de;
	chargedFragments[n].E           = e;
	chargedFragments[n].reconEnergy = recE;
	chargedFragments[n].x           = _x;
	chargedFragments[n].y           = _y;
}

// Polar angle of neutron fragment
void RootOutput::SetThetaNeut(double nth) {
	thetaNeut = nth;
	hist_neut_theta->Fill(nth);
}

// Secondary elastic scattering polar angle
// Assumes primary elastic scattering polar angle has already been set
void RootOutput::SetThetaElastS(double thEl) {
	thetaElastS = thEl;
	hist_theta_beam_S_sharp->Fill(sampler.thetaElastic);
  hist_theta_beam_S_recon->Fill(thEl);
}

// Primary relative energy (decay energy)
void RootOutput::SetErelP(double e) {
	ErelP = e;	
	hist_Erel_P->Fill(e);
}

// Secondary relative energy (decay energy)
void RootOutput::SetErelS(double e) {
	ErelS = e;
	hist_Erel->Fill(e);
}

// Excitation energy of fragment
void RootOutput::SetEx(double e) {
	Ex = e;
  hist_Ex->Fill(e);
}

// Cos(theta) of heavy decay fragment
// -1 or 1 is longitudinal decay, 0 is transverse decay
void RootOutput::SetCosThetaH(double c) {
	cosThetaH = c;
  cos_thetaH->Fill(c);
}

// Denotes if the beam fragment hit a detector or not
// Boolean values, so 1 is a hit and 0 is not
void RootOutput::SetIsElasticHit(int h) {
	isElasticHit = h;
}

// Denotes a successful detection of the decay fragments
// in the Si detectors
void RootOutput::SetIsFragDet(bool h) {
	isFragDet = h;
}

// Set values sampled from (in)elastic angular distributions
void RootOutput::SetSampledValues(SampledValues* s) {
	sampler = *s;

	hist_phi_P->Fill(sampler.phi);
	hist_theta_beam_P->Fill(sampler.thetaElastic);
	hist_theta_P->Fill(sampler.thetaLab);
	hist_vel_P->Fill(sampler.Vpplab);
	kinematic_circle->Fill(sampler.VppX, sampler.VppZ);
}




