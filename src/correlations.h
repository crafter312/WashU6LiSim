#ifndef _correlations
#define _correlations

#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>

#include "constants.h"
#include "random.h"

using namespace std;

struct SampledValues {
	// Elastic values
	double phi{-1};          // azimuthal angle sampled from uniform distribution
	double thetaElastic{-1}; // sampled elastic polar angle in lab frame

	// Inelastic values
	double thetaLab{-1};     // lab angle of outgoing projectile
	double Vpplab{-1};       // lab velocity of outgoing projectile
	double VppX{-1};         // lab velocity x component
	double VppY{-1};         // lab velocity y component
	double VppZ{-1};         // lab velocity z component

	void CalculateCartesian() {
		VppX = Vpplab * sin(thetaLab) * cos(phi); // x
		VppY = Vpplab * sin(thetaLab) * sin(phi); // y
		VppZ = Vpplab * cos(thetaLab);            // z
	};

	void Clear() {
		phi = -1;
		thetaElastic = -1;
		thetaLab = -1;
		Vpplab = -1;
		VppX = -1;
		VppY = -1;
		VppZ = -1;
	};
};

class Correlations {
  public:
    Correlations(string*, string, double, double, double*, double*, size_t);
    ~Correlations();

    CRandom ran;
    
    // Kinematic values
		double E;            // incoming beam energy
    double Exp;          // projectile excitation energy
		double Ext;          // target excitation energy
    double thetaCM;      // sampled inelastic polar angle in CM frame
    double thetaTarg;    // lab angle of outgoing target
    double Vttlab;       // lab velocity of outgoing target

		// Kinematic values to be saved
		SampledValues sampledValues;

		// Inelastic exit channel input
		string *filenames;  // inelastic input file names
		int *lengths;       // number of data points for each input file
		double **th_file;   // angular array for each input file
    double **Xsec_file; // differential cross section array for each input file
    
		// Elastic exit channel input
		string fileelastic;
    int lenElastic;
    double *th_elastic;
    double *Xsec_elastic;

		int nexits;    // total number of inelastic exit channels
		double *Xsecs; // total exit channel cross sections
		double *Exts;  // target excitation for each exit channel

    void randomAngles();             // samples elastic and inelastic scattering angles from input distributions
		void readelastic();              // reads elastic differential cross section Fresco file
    void readinelastic(string, int); // reads inelastic differential cross section Fresco file for given exit channel

	// Define constants used for calculations here
	private:
		double Mp;       // mass of projectile
		double Mt;       // mass of target
		double Mpp;      // mass of outgoing projectile (exit channel)
		double Mtt;      // mass of outgoing target (exit channel)
		double Mred;     // reduced mass
		double EnergyPA; // energy per nucleon
		double Vbeam;    // beam velocity
		double VCM;      // CM velocity
		double VpCM;     // velocity of projectile in CM frame
		double VtCM;     // velocity of target in CM frame
		double ECMin;    // kinetic energy of incoming target and projectile in CM frame
		double Qrxn;     // Q value of projectile decay

		void setConstants();       // calculates angle-independent kinematic values
		void calculateLabAngles(); // calculates lab frame velocities and angles for each sampled inelastic angle in CM frame
		
};

#endif
