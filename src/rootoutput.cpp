/**
 * Author: Henry Webb (h.s.webb@wustl.edu)
 * Created: 22 May 2024
 */

#include "rootoutput.h"

#include "constants.h"

// Input:
//	nFrags -- number of decay fragments
RootOutput::RootOutput(string suffix, int n, bool hasNeutron) : nFrags(n - hasNeutron) {
	realFragments.resize(nFrags);
	reconFragments.resize(nFrags);
	//for (int i = 0; i < nFrags; i++)
	//	chargedFragments[i] = CFragS();

	// Define output file
  string fileout = "sim" + suffix + ".root";
  file = new TFile(fileout.c_str(), "RECREATE");
	file->cd();

	/**** Initialize tree ****/

	t = new TTree("t", "t");
	t->Branch("realFragments", &realFragments);
	t->Branch("reconFragments", &reconFragments);
	t->Branch("reconElastic", &elastic);
	t->Branch("ErelPSampled", &ErelP);
	t->Branch("ErelPReal", &ErelPRecon);
	t->Branch("Ex", &Ex);
	t->Branch("cosThetaH", &cosThetaH);
	t->Branch("isElasticHit", &isElasticHit);
	t->Branch("isFragDet", &isFragDet);

	// Conditional neutron branches
	if (hasNeutron) {
		t->Branch("isNeutHit", &isNeutHit);
		t->Branch("neutron", &neutron);
	}

	// Reworked simulation branches
	t->Branch("sampler", &sampler);
	t->Branch("recon", &recon);

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
	// Fill histograms from tree
	t->Draw("recon.energy>>Erel","!std::isnan(recon.energy)","goff");
	t->Draw("recon.velocity>>vel_S","!std::isnan(recon.velocity)","goff");
	t->Draw("recon.theta>>theta_S","!std::isnan(recon.theta)","goff");
	t->Draw("recon.phi>>phi_S","!std::isnan(recon.phi)","goff");

	// Cleanup
	file->Write();
	file->Close();
}

void RootOutput::Fill() {
	t->Fill();
}

void RootOutput::Clear() {
	for (int i = 0; i < nFrags; i++) {
		realFragments[i].clear();
		reconFragments[i].clear();
	}
	elastic.clear();
	ErelP = NAN;
	ErelPRecon = NAN;
	Ex = NAN;
	cosThetaH = NAN;
	isElasticHit = 0;
	isFragDet = false;

	sampler.Clear();
	recon.Clear();

	isNeutHit = false;
	neutron.clear();
}

// Input:
//	n -- the fragment index (not including neutron)
//	de -- the energy lost in the thin front Si detector
//	e -- the remaining energy lost in the thick back Si detector
//	recE -- reconstructed fragment energy after detection
//	_x -- x position of fragment in detector
//	_y -- y position of fragment in detector
void RootOutput::SetRealFragment(int n, double de, double e, double recE, double _x, double _y, double _theta) {
	realFragments[n].DE          = de;
	realFragments[n].E           = e;
	realFragments[n].reconEnergy = recE;
	realFragments[n].x           = _x;
	realFragments[n].y           = _y;
	realFragments[n].thetaLab    = _theta;
}

// Input:
//	n -- the fragment index (not including neutron)
//	de -- the energy lost in the thin front Si detector
//	e -- the remaining energy lost in the thick back Si detector
//	recE -- reconstructed fragment energy after detection
//	_x -- x position of fragment in detector
//	_y -- y position of fragment in detector
void RootOutput::SetReconFragment(int n, double de, double e, double recE, double _x, double _y, double _theta) {
	reconFragments[n].DE          = de;
	reconFragments[n].E           = e;
	reconFragments[n].reconEnergy = recE;
	reconFragments[n].x           = _x;
	reconFragments[n].y           = _y;
	reconFragments[n].thetaLab    = _theta;
}

// Input:
//	de -- the energy lost in the thin front Si detector
//	e -- the remaining energy lost in the thick back Si detector
//	recE -- reconstructed fragment energy after detection
//	_x -- x position of fragment in detector
//	_y -- y position of fragment in detector
void RootOutput::SetElastic(double de, double e, double recE, double _x, double _y, double _theta) {
	elastic.DE          = de;
	elastic.E           = e;
	elastic.reconEnergy = recE;
	elastic.x           = _x;
	elastic.y           = _y;
	elastic.thetaLab    = _theta;

	hist_theta_beam_S_sharp->Fill(sampler.thetaElastic);
  hist_theta_beam_S_recon->Fill(_theta);
}

// Time, energy, angle, and position of neutron fragment
void RootOutput::SetNeut(double t, double E, double th, double cos, double x, double y, double z) {
	neutron.t = t;
	neutron.E = E;
	neutron.thetaLab = th;
	neutron.cos_thetaN = cos;
	neutron.pos[0] = x;
	neutron.pos[1] = y;
	neutron.pos[2] = z;
	hist_neut_theta->Fill(th);
}

// Primary relative energy (decay energy)
void RootOutput::SetErelP(double e) {
	ErelP = e;	
	hist_Erel_P->Fill(e);
}

// Primary relative energy reconstructed from real fragments
void RootOutput::SetErelPRecon(double e) {
	ErelPRecon = e;
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
void RootOutput::SetIsElasticHit(bool h) {
	isElasticHit = h;
}

// Denotes a successful detection of the charged decay
// fragments in the Si detectors
void RootOutput::SetIsFragDet(bool h) {
	isFragDet = h;
}

// Denotes a successful detection of the neutron decay
// fragment in the TexNeut detector
void RootOutput::SetIsNeutDet(bool h) {
	isNeutHit = h;
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

// Set reconstructed values from CDecay class
void RootOutput::SetReconValues(KinematicValues* r) {
	recon = *r;
}



