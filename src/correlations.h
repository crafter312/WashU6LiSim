#ifndef _correlations
#define _correlations

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>

#include "constants.h"
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

		double phi;          // azimuthal angle sampled from uniform distribution in degrees

		// Elastic values
		double thetaElastic; // sampled elastic polar angle in lab frame in degrees

		// Inelastic values
		double Ext;          // target excitation energy
		double thetaLab;     // lab angle of outgoing projectile in degrees
		double Vpplab;       // total lab velocity of outgoing projectile
		double VppX;         // lab velocity x component
		double VppY;         // lab velocity y component
		double VppZ;         // lab velocity z component
};

class Correlations {
  public:
    Correlations(string*, string, double, double, double*, double*, size_t);
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
