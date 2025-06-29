#include <iostream>
#include <cmath>
#include "rootoutput.h"
#include "Li6sim_alphapn.h"

using namespace std;

int main(int argc, char *argv[]) {
	
	/**** INPUT ARGUMENTS ****/

	// Total incoming beam energy in MeV, also used for the Fresco simulation.
	// If this changes, make sure to redo the Fresco simulations!
	double Ebeam = 63;

	double Ex     = 5.366; // excitation energy of parent fragment in MeV
	double gamma  = 0.541; // width of excited state of parent fragment in MeV

	// Default physical experiment parameters
	double distanceFromTarget = 150;        // distance of Gobbi from the target in mm
	string suffix             = "alphapn"; // output file suffix

	// Check for command line arguments, set default values if none are given
	//	- arg. 1 = beam energy in MeV
	//	- arg. 2 = Gobbi distance from target in mm
	//	- arg. 3 = intrinsic state width in MeV
	if (argc >= 3) {
		Ebeam = stod(argv[1]);
		distanceFromTarget = stod(argv[2]);
	}
	else cout << "WARNING: DEFAULT INPUT PARAMETERS BEING USED" << endl;

	// Optional third argument for state width, can supply the first two
	// without this one if desired.
	if (argc == 4) gamma = stod(argv[3]);

	// Initialize main simulation class
	Li6sim_alphapn sim(Ebeam, distanceFromTarget, Ex, gamma, suffix);
	//sim.AddExtraSuffix("perfTarg_noResolution2");

	// See Li6sim.h for default experiment parameters, which can
	// be changed via "Set..." commands as desired here.

	// Complete initialization of simulation class
	sim.Init();
	
	// Doesn't print absolutely everything, but can be changed if you wish
	sim.PrintSettings();

	// Initialize output manager
	RootOutput output(sim.GetSuffix(), sim.GetNFrags(), true);

	/**** MAIN EVENT LOOP ****/

	int Nevents = 100000; // number of events to simulate
	for (int index=0;index<Nevents;index++) {
	
		// Progress updates
		if (index%1000 == 0) cout << '\xd' <<	index << " of " << Nevents << flush;

		// Execute contents of event loop
		sim.DoSingleEventPreNeutron(output);
		sim.DoSingleEventPostNeutron(output);
	}

	// Output final statistics
	sim.DoFinalThings(Nevents);

	// Smart pointers now used by simulation class, so no manual cleanup is
	// required. I left the following comments in, though, because they're funny.

	//clean up, clean up
  //everybody everywhere
  //clean up, clean up
	//everybody do your share

	//beep at me when finished (sadly doesn't work anymore)
	//cout << "\a";
}
