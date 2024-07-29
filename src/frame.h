// Modified by Henry Webb on 24 July 2024

#ifndef frame_
#define frame_

#include <math.h>

#include "constants.h"

/**
 * Simple class to store important kinematic values for easy
 * ROOT output to a TTree object and do some basic conversions
 */
class KinematicValues {
public:
	KinematicValues();
	~KinematicValues();

	void Clear();

	void Sph2CartV();        // converts V from spherical to cartesian coordinates
	void Sph2CartPC();       // converts PC from spherical to cartesian coordinates
	void Cart2Sph();         // calculates spherical angles from V components
	void CalcVMag();         // calculates total V from components
	void CalcPCMag();        // calculates total PC from components
	void CalcCartV();        // calculates V components from total PC, total V, and PC components
	void CalcCartPC();       // calculates PC components from total V, total PC, and V components
	void AddVVec(double*);   // adds supplied vector to V vector
	void AddPCVec(double*);  // adds supplied vector to PC vector
	double GetV2();          // returns V squared
	double GetPC2();         // returns PC squared
	double GetVDot(double*); // returns dot product of V vector with supplied vector

	double energy {NAN};   // fragment energy in MeV
	double velocity {NAN}; // fragment velocity in cm/ns
	double pcTot {NAN};    // total momentum*c of fragment
	double v[3];           // velocity vector of fragment
	double pc[3];          // momentum*c vector of fragment MeV
	double theta {NAN};    // polar angle of fragment in radians
	double phi {NAN};      // azimuthal angle of fragment in radians
	double x {NAN};
	double y {NAN};
};

/**
 * Short class to handle calculation of Relativistic and
 * Non-Relativistic kinematic values
 *
 * give velocity, energies, transforms to new reference frame
 * using either non-relativistic or Relativistics equations
 */
class CFrame {
public:
	CFrame(double);

	double getVelocity(bool*);
	double getVelocityNewton();
	double getVelocityRel();
	double getEnergy(bool*);
	//double getEnergyFromMomentum();
	double getEnergyNewton();
	double getEnergyRel();
	void transformVelocity(double*, bool*);
	void transformVelocityNewton(double*);
	void transformVelocityRel(double*);
	void getVelocityFromMom(bool*);
	void getVelocityFromMomNewton();
	void getVelocityFromMomRel();
	void getMomFromVelocity();

	// Pass-throughs for KinematicValue functions
	void AddVVec(double, double, double);
	void AddPCVec(double, double, double);
	void ScaleVVec(double);
	void CalcVMag();
	void CalcPCMag();
	void Cart2Sph();
	void CalcCartV();
	void Sph2CartV();
	double GetV2();

	// Setters and Getters
	void SetEnergy(double);
	void SetVelocity(double);
	void SetVelocityComps(double, double, double);
	void SetMomComps(double, double, double);
	void SetTheta(double);
	void SetPhi(double);
	void SetXY(double, double);
	double GetEnergy(); // not to be confused with getEnergy(bool*) from above
	double GetVelocity();
	double GetPC();
	double GetVComp(int);
	double GetPCComp(int);
	double GetTheta();
	double GetPhi();
	double GetX();
	double GetY();

	KinematicValues* GetKinematicValues();

	double A;         // mass number
	double mass;      // rest mass of fragment in MeV
	double totEnergy; // rest mass + kinetic energy

private:
	KinematicValues kinematicValues; // stores all kinematic values useful for simulation output

	double mass2; // mass^2 for calculation purposes
	double gamma; //!< relativistic gamma factor
};


#endif
