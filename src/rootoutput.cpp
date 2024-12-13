/**
 * Author: Henry Webb (h.s.webb@wustl.edu)
 * Created: 22 May 2024
 */

#include "rootoutput.h"

// Input:
//	nFrags -- number of decay fragments
RootOutput::RootOutput(string suffix, int n) : nFrags(n) {
	chargedFragments.resize(nFrags);
	//for (int i = 0; i < nFrags; i++)
	//	chargedFragments[i] = CFragS();

	// Define output file
  string fileout = "sim" + suffix + ".root";
  file = new TFile(fileout.c_str(), "RECREATE");
	file->cd();

	/**** Initialize tree ****/

	t = new TTree("t", "t");
	t->Branch("chargedFragments", &chargedFragments);
	t->Branch("ENeut", &ENeut);
  t->Branch("thetaNeut", &thetaNeut);
	t->Branch("ErelP", &ErelP);
	t->Branch("Ex", &Ex);
	t->Branch("cosThetaH", &cosThetaH);
	t->Branch("isElasticHit", &isElasticHit);
	t->Branch("isFragDet", &isFragDet);
	t->Branch("Ex_TEST", &Ex_TEST);
	t->Branch("Erel_TEST", &Erel_TEST);
	t->Branch("Ex_TEST_CM", &Ex_TEST_CM);
	t->Branch("Erel_TEST_CM", &Erel_TEST_CM);

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
	for (int i = 0; i < nFrags; i++)
		chargedFragments[i].clear();
	ENeut = NAN;
	thetaNeut = NAN;
	thetaElastS = NAN;
	ErelP = NAN;
	Ex = NAN;
	cosThetaH = NAN;
	isElasticHit = 0;
	isFragDet = false;

	Ex_TEST = NAN;
	Erel_TEST = NAN;
	Ex_TEST_CM = NAN;
	Erel_TEST_CM = NAN;

	sampler.Clear();
	recon.Clear();
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

// Energy of neutron fragment
void RootOutput::SetENeut(double E) {
	ENeut = E;
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

// Set reconstructed values from CDecay class
void RootOutput::SetReconValues(KinematicValues* r) {
	recon = *r;
}

// Set test values in output tree. This function should
// be modified to include any desired test output values.
void RootOutput::SetTestValues(double Ex, double Erel, double ExCM, double ErelCM) {
	Ex_TEST = Ex;
	Erel_TEST = Erel;
	Ex_TEST_CM = ExCM;
	Erel_TEST_CM = ErelCM;
}



