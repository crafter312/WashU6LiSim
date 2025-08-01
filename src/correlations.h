#ifndef _correlations
#define _correlations

#include <exception>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "frame.h"
#include "loss.h"
#include "random.h"

using namespace std;

class SampledValues {
	public:
		SampledValues();
		~SampledValues();

		void Clear();

		void CalculateCartesian();
		
		double GetPhiRad();
		double GetThetaElasticRad();
		double GetThetaLabRad();

		double phi {NAN};          // azimuthal angle sampled from uniform distribution in degrees

		// Elastic values
		double thetaElastic {NAN}; // sampled elastic polar angle in lab frame in degrees

		// Inelastic values
		double Ext {NAN};          // target excitation energy
		double thetaLab {NAN};     // lab angle of outgoing projectile in degrees
		double Vpplab {NAN};       // total lab velocity of outgoing projectile
		double VppX {NAN};         // lab velocity x component
		double VppY {NAN};         // lab velocity y component
		double VppZ {NAN};         // lab velocity z component
};

class Correlations {
  public:
    Correlations(vector<string>, string, double, double, vector<double>, size_t, string, float);
    ~Correlations();

    CRandom ran;
    
    // Kinematic values
		double E;            // incoming beam energy
    double Exp;          // projectile excitation energy
    double thetaCM;      // sampled inelastic polar angle in CM frame
    double thetaTarg;    // lab angle of outgoing target
    double Vttlab;       // lab velocity of outgoing target

		// Kinematic values to be saved
		SampledValues sampledValues;

		// Inelastic exit channel input
		vector<string> filenames;  // inelastic input file names
		int *lengths;       // number of data points for #include "loss.h"each input file
		double **th_file;   // angular array for each input file
    double **Xsec_file; // differential cross section array for each input file
    
		// Elastic exit channel input
		string fileelastic;
    int lenElastic;
    double *th_elastic;
    double *Xsec_elastic;

		int nexits;           // total number of inelastic exit channels
		double *Xsecs;        // total exit channel cross sections
		vector<double> Exts;  // target excitation for each exit channel

    void randomAngles(double);         // samples elastic and inelastic scattering angles from input distributions
		void readelastic(float);           // reads elastic differential cross section Fresco file
    double readinelastic(string, int); // reads inelastic differential cross section Fresco file for given exit channel, returns integrated cross section

	// Define constants used for calculations here
	private:
		double Mp;             // mass of projectile (amu)
		double Mt;             // mass of target (amu)
		double Mpp;            // mass of outgoing projectile (exit channel) (amu)
		double Mtt;            // mass of outgoing target (exit channel) (amu)
		double mpp;            // mass of outgoing projectile (exit channel) (MeV / c^2)
		double mtt;            // mass of outgoing target (exit channel) (MeV / c^2)
		double mpptt;          // total mass of outgoing projectile and target (exit channel) (MeV / c^2)
		double EnergyPostLoss; // beam energy after energy loss in target
		double EnergyPA;       // energy per nucleon
		double gammaVbeam;     // relativistic beam momentum divided by mass
		double Vbeam;          // beam velocity
		double VCM;            // CM velocity
		double VCMvec[3] = { 0., 0., 0. };
		double VpCM;           // velocity of projectile in CM frame
		double VtCM;           // velocity of target in CM frame
		double ECMin;          // kinetic energy of incoming target and projectile in CM frame
		double Qrxn;           // Q value of projectile decay

		CLoss* ploss_C;        // loss file for incoming projectile in C target
		CFrame* framep;        // incoming beam frame
		CFrame* framet;        // incoming target frame
		CFrame* framepp;       // outgoing parent frame
		CFrame* framett;       // outgoing target frame

		void setConstants();             // calculates angle-independent kinematic values
		void calculateLabAngles(double); // calculates lab frame velocities and angles for each sampled inelastic angle in CM frame
		
};

#endif
