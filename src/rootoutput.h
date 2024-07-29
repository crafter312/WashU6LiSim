/**
 * Author: Henry Webb (h.s.webb@wustl.edu)
 * Created: 22 May 2024
 */

#ifndef _rootoutput
#define _rootoutput

#include <math.h>
#include <optional>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2S.h"

#include "constants.h"
#include "correlations.h"
#include "frame.h"

using namespace std;

/**
 * Lab velocity, angles, etc. of the parent fragment.
 * This can either be the primary distributions after sampling
 * or the secondary reconstructed distributions.
 */
struct PFragS {
	double vel { NAN };
	double velx { NAN };
	double vely { NAN };
	double velz { NAN };
	double phi { -1 };
	double theta { -1 };
	double Erel { -1 };

	void clear() {
		vel = NAN;
		velx = NAN;
		vely = NAN;
		velz = NAN;
		phi = -1;
		theta = -1;
		Erel = -1;
	};
};

/**
 * Energy, position, etc. of charged fragment.
 * This is either alpha or proton in my case.
 */
struct CFragS {
	double DE { -1 };
	double E { -1 };
	double reconEnergy { -1 };
	double x { NAN };
	double y { NAN };

	void clear() {
		DE = -1;
		E = -1;
		reconEnergy = -1;
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

		void SetSecondary(double, double, double);
		void SetFragment(int, double, double, double, double, double);
		void SetThetaNeut(double);
		void SetThetaElastS(double);
		void SetErelP(double);
		void SetErelS(double);
		void SetEx(double);
		void SetCosThetaH(double);
		void SetIsElasticHit(int);
		void SetIsFragDet(bool);

		void SetSampledValues(SampledValues*);

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

		PFragS parentSecondary {};
		CFragS* chargedFragments;
		int nFrags;
		double thetaNeut { -1 };
		double thetaElastS { -1 };
		double ErelP { -1 };
		double ErelS { -1 };
		double Ex { -1 };
		double cosThetaH { NAN };
		int isElasticHit { -1 };
		bool isFragDet { 0 };

		SampledValues sampler {};
		
};

#endif
