#ifndef _Li6sim_alphapn
#define _Li6sim_alphapn

#include "correlations.h"
#include "decay.h"
#include "frag.h"
#include "Gobbiarray.h"
#include "rootoutput.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

class Li6sim_alphapn {
public:
	Li6sim_alphapn(double, double, double, double, std::string);
	Li6sim_alphapn(double, double, double, double, std::string, float, float);
	~Li6sim_alphapn();

	// Functions to set simulation-specific variables
	void SetTargetThickness(float th) { thickness = th; }
	void SetNeutTRes(float res) { neutTRes = res; }
	void SetGobbiRes(float res) { GobbiRes = res; }
	void SetBeamSpotRadius(float r) { targetSize = r; }
	void SetSmallAngleScatteringScale(float s) { scale = s; }
	void SetRelativistic(bool b) { einstein = b; }
	void SetUseRealP(bool b) { useRealP = b; }
	void SetEnableMARS(bool b) { mars = b; }

	size_t GetNFrags() const { return Nfrag; }
	std::string GetSuffix() const { return suffix; }

	// Extra functions
	void AddExtraSuffix(std::string s) { suffix += "_" + s; }

	// Simulation functions (these split what previously was the main simulation file
	// into multiple functions, for ease of implementing into an external program)

	// This performs all remaining initial routines, and returns a non-empty string
	// in the case where something has gone wrong. The responsibility is on the
	// function caller to make sure that the program exits correctly if a non-empty
	// string is returned.
	std::string Init();

	// This handles all the cout statements of the simulation parameters
	void PrintSettings();

	// Together, these functions handle the contents of the event loop, taking a
	// reference to RootOutput. Again, a non-empty string is returned if something
	// has gone wrong, and it is the responsibility of the function caller to make
	// sure the program exits properly if this happens.

	// This function handles all of the elements of the simulation which do not
	// make use of the neutron results from the decay calculation.
	std::string DoSingleEventPreNeutron(RootOutput&);

	// This function handles all simulation elements that require neutron results,
	// which can be set from an external simulation if desired.
	std::string DoSingleEventPostNeutron(RootOutput&);

	// Handle final things like outputting statistics and such
	void DoFinalThings(int);

private:

	// Simulation-specific variables
	double Ebeam{-1.};
	double distanceFromTarget{-1.};
	double Ex{-1.};
	double gamma{-1.};
	std::string suffix{""};
	double Q{NAN};

	// Simulation-specific variables with default values
	float thickness{3.026};     // target thickness in mg/cm^2 (copied from Nic's experiment)
	float neutTRes{0.5};        // neutron timing resolution (sigma) in ns
	float GobbiRes{0.02};       // Si-Si resolution;
	float targetSize{1.};       // diameter of beam spot size in mm
	float scale{1.38};          // scales the magnitude of small angle scattering
	bool einstein{true};        // switch for newtonian(false) or relativistic(true) kinematics
	bool mars{false};           // switch to enable MARS momentum acceptance

	bool useRealP{false}; // true means use real angle and energies of fragment
	                      // for event reconstruction, to check effect of
	                      // detector resolution

	// These two can be changed using one of the constructors
	float b{8.0};               // mm beam axis to Gobbi frame dimension,
	float RadiusCollimator{0.}; // mm Gobbi collimator outer radius

	// Files for energy loss in C target
	std::string Loss_Li_in_C{"Lithium_C.loss"};
	std::string Loss_He_in_C{"Helium_C.loss"};
	std::string Loss_p_in_C{"Hydrogen_C.loss"};
	std::string Loss_n_in_C{"Hydrogen_C.loss"}; // LOSS FILE NOT USED, TEMPORARY STAND-IN

	// Files for energy loss in Si detector
	std::string Loss_Li_in_Si{"Lithium_Si.loss"};
	std::string Loss_He_in_Si{"Helium_Si.loss"};
	std::string Loss_p_in_Si{"Hydrogen_Si.loss"};
	std::string Loss_n_in_Si{"Hydrogen_Si.loss"}; // LOSS FILE NOT USED, TEMPORARY STAND-IN

	// Production reaction exit channel variables
	size_t nexits;                      // number of target excited state exit channels for production reaction
	std::vector<double> Exts;           // outgoing target excitation energy for each exit channel
	std::vector<std::string> Xsecfiles; // file paths for differential cross-sections of each exit channel
	std::string elasXsecfile;           // file path for differential cross-section of elastic scattering

	// Decay fragments. Use smart pointers so that
	// objects can be initialized at a later point.
	size_t Nfrag;
	std::vector<std::shared_ptr<CFrag>> frag;
	std::unique_ptr<CFrag> fragBeam;

	// Simulation classes
	std::shared_ptr<Gobbiarray> gobbi;
	std::unique_ptr<CDecay> decay;
	std::unique_ptr<Correlations> sampler;

	// Physics variables
	double pc0{-1.};         // pc of beam
	double P_acceptance{0.}; // momentum acceptance of MARS, if enabled

	// Loop counter variables
	int Npunch{0};
	int Nmiss{0};
	int N2Hit{0};
	int Nstuck{0};
	int Ndet{0};
	int Nbeamscat{0};
};

#endif
