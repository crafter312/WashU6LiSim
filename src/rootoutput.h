/**
 * Author: Henry Webb (h.s.webb@wustl.edu)
 * Created: 22 May 2024
 */

#ifndef _rootoutput
#define _rootoutput

#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2S.h"

#include "correlations.h"
#include "frame.h"

using namespace std;

/**
 * Energy, position, etc. of charged fragment.
 * This is either alpha or proton in my case.
 */
struct CFragS {
	double DE{NAN};
	double E{NAN};
	double reconEnergy{NAN};
	double x{NAN};
	double y{NAN};
	double thetaLab{NAN};

	void clear() {
		DE = NAN;
		E = NAN;
		reconEnergy = NAN;
		x = NAN;
		y = NAN;
		thetaLab = NAN;
	};
};

/**
 * Energy, position, etc. of neutron fragment.
 */
struct NeutFrag {
	bool wasDarkScattered{false};
	double t{NAN};
	double E{NAN};
	double thetaLab{NAN};
	double cos_thetaN{NAN};
	double pos[3]{NAN, NAN, NAN};

	void clear() {
		wasDarkScattered = false;
		t = NAN;
		E = NAN;
		thetaLab = NAN;
		cos_thetaN = NAN;
		pos[0] = NAN;
		pos[1] = NAN;
		pos[2] = NAN;
	};
};

/**
 * RootOutput class.
 * This handles the output file, histograms, and tree from the simulation.
 */
class RootOutput {
	public:
  	TH2I* DEE;
  	TH1F* protonenergy;
  	TH1F* alphaenergy;
  	TH2F* hist_Erel_thetaH;
  	TH1I* hist_Ex_trans;
  	TH1I* hist_Ex_trans_narrow;
		TH2I* hist_Ex_DE;
		TH2S* protonXY_S;
  	TH2S* coreXY_S;

		RootOutput(string, int, bool = true);
		~RootOutput();

		void Fill();
		void Clear();

		void SetRealFragment(int, double, double, double, double, double, double);
		void SetReconFragment(int, double, double, double, double, double, double);
		void SetElastic(double, double, double, double, double, double);
		void SetNeut(bool, double, double, double, double, double, double, double);
		void SetErelP(double);
		void SetErelPRecon(double);
		void SetEx(double);
		void SetCosThetaH(double);
		void SetIsElasticHit(bool);
		void SetIsFragDet(bool);
		void SetIsNeutDet(bool);
		void SetTargetEloss(double, double, double, double, double);

		void SetSampledValues(SampledValues*);
		void SetReconValues(KinematicValues*);
		void SetDETests(std::vector<double>&);

	private:
		TFile* file;
		TTree* t;

		TH1F* hist_vel_P;
		TH1F* hist_theta_P;
		TH1F* hist_phi_P;
  	TH2F* kinematic_circle;
		TH1F* hist_theta_beam_P;
		TH1I* hist_Erel_P;
		TH1F* hist_vel_S;
  	TH1F* hist_theta_S;
  	TH1F* hist_phi_S;
		TH1F* hist_theta_beam_S_sharp;
  	TH1F* hist_theta_beam_S_recon;
		TH1F* hist_neut_theta;
  	TH1I* hist_Erel;
  	TH1I* hist_Ex;
  	TH1F* cos_thetaH;

		std::vector<CFragS> realFragments;
		std::vector<CFragS> reconFragments;
		std::vector<double> dETests;
		CFragS elastic;
		int nFrags;
		double ErelP {NAN};         // this is the sampled Breit-Wigner distribution directly from the decay class
		double ErelPRecon {NAN};    // this is the decay energy calculated from the real fragments immediately after the decay
		double Ex {NAN};
		double cosThetaH {NAN};
		bool isElasticHit {false};  // flag for valid elastic scattering hit in Gobbi
		bool isFragDet {false};     // flag for valid alpha + p coincidence in Gobbi
		double targEloss {NAN};     // real total energy loss of beam and fragments in the target
		double targElossRec {NAN};  // total energy loss of beam and fragments in the target + diamond detector resolution
		double inthick {NAN};       // reaction position in target from upstream side, in mg/cm^2
		double inthickrec {NAN};    // reconstructed reaction position in target from upstream side, in mg/cm^2
		double inthickrecimp {NAN}; // reconstructed reaction position in target from upstream side using improved procedure, in mg/cm^2

		// Conditional neutron values
		bool isNeutHit {false}; // flag for valid neutron hit in TexNeut
		NeutFrag neutron;

		SampledValues sampler {};
		KinematicValues recon {};
		
};

#endif
