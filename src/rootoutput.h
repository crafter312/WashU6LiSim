/**
 * Author: Henry Webb (h.s.webb@wustl.edu)
 * Created: 22 May 2024
 */

#ifndef _rootoutput
#define _rootoutput

#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2S.h"

#include "constants.h"
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

	void clear() {
		DE = NAN;
		E = NAN;
		reconEnergy = NAN;
		x = NAN;
		y = NAN;
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

		RootOutput(string, int);
		~RootOutput();

		void Fill();
		void Clear();

		void SetFragment(int, double, double, double, double, double);
		void SetENeut(double);
		void SetThetaNeut(double);
		void SetThetaElastS(double);
		void SetErelP(double);
		void SetEx(double);
		void SetCosThetaH(double);
		void SetIsElasticHit(int);
		void SetIsFragDet(bool);

		void SetTestValues(double, double);

		void SetSampledValues(SampledValues*);
		void SetReconValues(KinematicValues*);

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

		std::vector<CFragS> chargedFragments;
		int nFrags;
		double ENeut {NAN};
		double thetaNeut {NAN};
		double thetaElastS {NAN};
		double ErelP {NAN};
		double Ex {NAN};
		double cosThetaH {NAN};
		int isElasticHit {0};
		bool isFragDet {false};

		// Test outputs
		double Ex_TEST {NAN};
		double Erel_TEST {NAN};

		SampledValues sampler {};
		KinematicValues recon {};
		
};

#endif
