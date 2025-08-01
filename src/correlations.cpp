// Addition added to sim.cpp to incorparate differential cross section into the
// selection of angle theta
// Author: Nicolas Dronchi, 07/23/2020
// Modified: Henry Webb, 05/15/2024 07/22/2024 11/19/2024

#include "correlations.h"

#include "constants.h"

#include <iostream>
#include <stdexcept>

using namespace std;

SampledValues::SampledValues() {}

SampledValues::~SampledValues() {}

void SampledValues::Clear() {
	phi = NAN;
	thetaElastic = NAN;
	thetaLab = NAN;
	Vpplab = NAN;
	VppX = NAN;
	VppY = NAN;
	VppZ = NAN;
}

void SampledValues::CalculateCartesian() {
	VppX = Vpplab * sin(thetaLab) * cos(phi); // x
	VppY = Vpplab * sin(thetaLab) * sin(phi); // y
	VppZ = Vpplab * cos(thetaLab);            // z
}

double SampledValues::GetPhiRad() {
	return phi * deg_to_rad;
}

double SampledValues::GetThetaElasticRad() {
	return thetaElastic * deg_to_rad;
}

double SampledValues::GetThetaLabRad() {
	return thetaLab * deg_to_rad;
}

/**********************************************************************************************/

Correlations::Correlations(vector<string> filenamein, string fileelasticin, double E0, double Ex0, vector<double> Ext0s, size_t n, string lossfile_C, float thickness) {
	if((filenamein.size() != n) || (Ext0s.size() != n))
		throw invalid_argument("vector filenamein and Ext0s lengns ");

	CRandom ran;
	filenames = filenamein;
	fileelastic = fileelasticin;

	nexits = n;

	E    = E0;	  // MeV
	Exp  = Ex0;   // MeV
	Exts = Ext0s; // MeV

	// Declare arrays for inelastic input
	lengths   = new int[nexits];
	th_file   = new double*[nexits];
	Xsec_file = new double*[nexits];

	// Some constants and calculations are same for all sampled angles
	setConstants();

	// Initialize loss object from file, must be done after setConstants
	ploss_C = new CLoss(lossfile_C, Mp);

	// Read the Fresco fort.201 file containing unscaled elastic cross 
	// section data in the center of mass
	readelastic(thickness);

	// Read all the supplied Fresco fort.202, fort.203, etc. files containing
	// inelastic cross section data in the center of mass. Also record
	// total integrated cross sections for sampling.
	Xsecs = new double[nexits];
	for (int i = 0; i < nexits; i++) {
		Xsecs[i] = readinelastic(filenames[i], i);
		if (i == 0) continue;
		Xsecs[i] += Xsecs[i - 1];
	}
	cout << "Total inelastic scattering cross section (mb): " << Xsecs[nexits - 1] << endl;

	// Probability distribution to pick an exit channel using the same logic
	// that's used to create a probability distribution from the differential cross 
	// sections. This is created using integrated total cross sections from
	// the input differential cross section files. The number and position of
	// these differential cross section files in the input must match the input
	// excitation energies!
	for (int i = 0; i < nexits; i++)
		Xsecs[i] /= Xsecs[nexits - 1];
}

Correlations::~Correlations() {
	delete []th_file;
	delete []Xsec_file;
	delete []th_elastic;
  delete []Xsec_elastic;

	delete framep;
	delete framet;
	delete framepp;
	delete framett;
}

void Correlations::randomAngles(double thick) {
	sampledValues.phi = 2. * pi * ran.Rndm();

	/**** INELASTIC ****/
	// Probability to accept a random theta is dependent on the cross section

	// First pick an exit channel
	float probtr = ran.Rndm();
	int i = 0;
	for (;;) {
		if ((Xsecs[i] > probtr) || (i + 1 == nexits)) break;
		i++;
	}
	sampledValues.Ext = Exts[i];
	double* Xsec = Xsec_file[i];
	double* th = th_file[i];
	int length = lengths[i];

	// Iterate through scaled/integrated cross section until random threshold is met
	probtr = ran.Rndm();
	int j = 0;
	for (;;) {
		if ((Xsec[j] > probtr) || (j + 1 == length)) break;
		j++;
	}

	// Pick random angle between chosen angle and next one
	thetaCM = (j + 1 == length) ? th[length - 1] : th[j] + (th[j + 1] - th[j]) * ran.Rndm();

	/**** ELASTIC ****/

	probtr = ran.Rndm();

	// Iterate through scaled/integrated cross section until random threshold is met
	int ii = 0;
	for (;;) {
		if ((Xsec_elastic[ii] > probtr) || (ii + 1 == lenElastic)) break;
		ii++;
	}
	sampledValues.thetaElastic = (ii + 1 == length) ? th_elastic[lenElastic - 1] : th_elastic[ii] + (th_elastic[ii + 1] - th_elastic[ii]) * ran.Rndm();
	
	calculateLabAngles(thick);
}

void Correlations::setConstants() {
	Mp     = Mass_7Li / m0;           // mass of projectile Li7 (amu)
	Mt     = Mass_12C / m0;           // mass of target C12 (amu)
	Mpp    = Mass_6Li / m0;           // mass of outgoing projectile Li6 (amu)
	Mtt    = Mass_13C / m0;           // mass of outgoing target C13 (amu)
	mpp    = Mass_6Li;                // mass of outgoing projectile Li6 (MeV / c^2)
	mtt    = Mass_13C;                // mass of outgoing target C13 (MeV / c^2)
	mpptt  = mpp + mtt;               // total mass of outgoing projectile and target (exit channel) (MeV / c^2)

	
	framep = new CFrame(Mp); // Incoming beam frame
	framet = new CFrame(Mt); // Incoming target frame
	framepp = new CFrame(Mpp); // Outgoing parent frame
	framett = new CFrame(Mtt); // Outgoing target frame

	Qrxn = mass_7Li + mass_12C - mass_6Li - mass_13C; // Q value for 7Li + 12C -> 6Li + 13C, ~ -2.30478473 MeV
	cout << "n-transfer Q-value: " << Qrxn << endl;
}

void Correlations::calculateLabAngles(double thick) {
	// Reset frame variables
	framep->SetTheta(0.);
	framep->SetPhi(0.);
	framet->SetTheta(0.);
	framet->SetPhi(0.);
	framet->SetEnergy(0.);
	framet->getVelocityRel();

	// Account for energy loss in target
	EnergyPostLoss = ploss_C->getEout(E, thick);
	framep->SetEnergy(EnergyPostLoss);
	framep->getVelocityRel();

	// Calculate CM velocity
	VCM = framep->GetPC() * c / (framep->totEnergy + framet->totEnergy);
	VCMvec[2] = -VCM;
	
	// Then transform to CM frame
	framet->transformVelocityRel(VCMvec);
	framep->transformVelocityRel(VCMvec);

	ECMin = framet->GetEnergy() + framep->GetEnergy(); // kinetic energy of incoming target and projectile in CM frame

	// Energy of projectile and target in exit channel
	// Note that the target nucleus excitation energy 
	// changes depending on the exit channel chosen
	double ECMout = ECMin + Qrxn - Exp - sampledValues.Ext;
	double ECMout2 = ECMout*ECMout;

	// Resulting momentum of both fragments in exit channel
	// This comes from setting up a standard relativistic
	// equation for total kinetic energy and solving for pc
	// using Wolfram Alpha
	double PCCMout = 0.5 * sqrt((2. * mpp * ECMout) + ECMout2) * sqrt(((2. * mtt) + ECMout) * ((2. * mpptt) + ECMout)) / (mpptt + ECMout);

	// Set values of parent fragment in CM frame
	framepp->SetTheta(thetaCM);
	framepp->SetPhi(sampledValues.phi);
	framepp->totEnergy = sqrt((mpp*mpp) + (PCCMout*PCCMout));
	framepp->SetVelocity(PCCMout * c / framepp->totEnergy);
	framepp->Sph2CartV();

	// Transform outgoing parent fragment to lab frame
	VCMvec[2] *= -1; // change sign of velocity to boost instead of "unboost"
	framepp->transformVelocityRel(VCMvec);

	/**** lab velocity and angle of the projectile ****/
	sampledValues.Vpplab = framepp->GetVelocity();
	sampledValues.thetaLab = framepp->GetTheta();
	sampledValues.CalculateCartesian();

	// Set values of outgoing target in CM frame (180 deg from frag)
	framett->SetTheta(pi - thetaCM);
	framett->SetPhi(sampledValues.phi - (pi * (1 - (2 * (sampledValues.phi < pi))))); // -pi for phi>=pi, +pi for phi<pi
	framett->totEnergy = sqrt((mtt*mtt) + (PCCMout*PCCMout));
	framett->SetVelocity(PCCMout * c / framett->totEnergy);
	framett->Sph2CartV();

	// Transform outgoing parent fragment to lab frame
	framepp->transformVelocityRel(VCMvec);

	/**** lab velocity and angle of the target ****/
	Vttlab = framett->GetVelocity();
	thetaTarg = framett->GetTheta();

	// Convert angles to degrees
	sampledValues.phi *= rad_to_deg;
	sampledValues.thetaElastic *= rad_to_deg;
	sampledValues.thetaLab *= rad_to_deg;
}

// Elastic scattering angles are converted to the lab frame when reading the file
void Correlations::readelastic(float thickness) {
	// Open file
	fstream file;
	cout << "Elastic Differential Cross Section filename " << fileelastic << endl;
	file.open(fileelastic, ios::in);
	if (!file.is_open())
		throw invalid_argument("Error opening file " + fileelastic);

	string line;

	// Extract # angles from first line
	getline(file, line);
	lenElastic = stoi(line.substr(1, line.find(" ")));
	cout << "# elastic angles: " << lenElastic << endl;

	// Declare data arrays
	th_elastic = new double[lenElastic];
	Xsec_elastic = new double[lenElastic];

	// Reset frame variables
	framep->SetTheta(0.);
	framep->SetPhi(0.);
	framet->SetTheta(0.);
	framet->SetPhi(0.);
	framet->SetEnergy(0.);
	framet->getVelocityRel();

	// Account for energy loss in target
	EnergyPostLoss = ploss_C->getEout(E, thickness*0.5);
	framep->SetEnergy(EnergyPostLoss);
	framep->getVelocityRel();

	// Calculate CM velocity
	VCM = framep->GetPC() * c / (framep->totEnergy + framet->totEnergy);
	VCMvec[2] = -VCM;
	
	// Then transform to CM frame
	framet->transformVelocityRel(VCMvec);
	framep->transformVelocityRel(VCMvec);

	ECMin = framet->GetEnergy() + framep->GetEnergy(); // kinetic energy of incoming target and projectile in CM frame
	double ECMout = ECMin;                             // energy out = energy in (elastic scattering)
	double ECMout2 = ECMout*ECMout;

	// Resulting momentum of both fragments in exit channel
	// This comes from setting up a standard relativistic
	// equation for total kinetic energy and solving for pc
	// using Wolfram Alpha
	double PCCMout = 0.5 * sqrt((2. * mpp * ECMout) + ECMout2) * sqrt(((2. * mtt) + ECMout) * ((2. * mpptt) + ECMout)) / (mpptt + ECMout);

	VCMvec[2] *= -1; // change sign of velocity to boost instead of "unboost"

	// Get data from file
	int count = 0;
	double theta_i, theta_rad, thetaLab;
	while (getline(file, line)) {
		if ((line.compare(0,1,"#") == 0) || (line.compare(0,1,"@") == 0) || (line.find("END") != string::npos))
			continue;

		stringstream sstr(line);
		sstr >> theta_i >> Xsec_elastic[count]; // extracts angle and cross section values from a particular input file line
		theta_rad = theta_i * deg_to_rad;
		Xsec_elastic[count] *= sin(theta_rad);

		/**** TRANSFORM ANGLE TO LAB FRAME ****/

		// Set values of parent fragment in CM frame
		framepp->SetTheta(theta_rad);
		framepp->SetPhi(0.);
		framepp->totEnergy = sqrt((mpp*mpp) + (PCCMout*PCCMout));
		framepp->SetVelocity(PCCMout * c / framepp->totEnergy);
		framepp->Sph2CartV();

		// Transform outgoing parent fragment to lab frame
		framepp->transformVelocityRel(VCMvec);
		thetaLab = framepp->GetTheta();
		th_elastic[count] = thetaLab;

		// Integrate differential cross section for theta > 3 deg.
		Xsec_elastic[count] = (thetaLab > (3. * deg_to_rad)) ? Xsec_elastic[count] + Xsec_elastic[count - 1] : 0;
		count++;
	}

	// Total elastic cross section for thetaLab > 3 degrees using left Reimann sum
	cout << "Elastic cross section (mb): " << twopi * Xsec_elastic[lenElastic - 2] * 180. / (double)(lenElastic - 1) * deg_to_rad << endl;

	// scale probability distribution to have range (0, 1)
	double norm = Xsec_elastic[lenElastic - 1];
	for (int i = 0; i < lenElastic; i++)
		Xsec_elastic[i] /= norm;

	file.close();
}

double Correlations::readinelastic(string filename, int ind) {
	fstream file;
	cout << "Inelastic Differential Cross Section file: " << filename << endl;
	file.open(filename, ios::in);
	if (!file.is_open())
		throw invalid_argument("Error opening file " + filename);

	string line;

	// Extract # angles from first line
	getline(file, line);
	lengths[ind] = stoi(line.substr(1, line.find(" ")));
	cout << "# inelastic angles: " << lengths[ind] << endl;

	// create arrays to store the cross section data
	th_file[ind] = new double[lengths[ind]];
	Xsec_file[ind] = new double[lengths[ind]];
	double* th = th_file[ind];
	double* Xsec = Xsec_file[ind];

	// go through file and take out cross section data
	int count = 0;
	double theta_i;
	while (getline(file, line)) {
		if ((line.compare(0,1,"#") == 0) || (line.compare(0,1,"@") == 0) || (line.find("END") != string::npos))
			continue;

		// store cross section data in arrays
		stringstream sstr(line);
		sstr >> theta_i >> Xsec[count];
		th[count] = theta_i*deg_to_rad;
		Xsec[count] *= sin(th[count]);

		// integrate differential cross section
		if (count == 0) {
			count++;
			continue;
		}

		Xsec[count] += Xsec[count-1];
		count++;
	}

	// Total inelastic cross section using left Reimann sum
	// The first bin is slightly narrower, so this is not totally correct
	double xsecTot = twopi * Xsec[lengths[ind] - 2] * 180. / (double)(lengths[ind] - 1) * deg_to_rad;
	cout << "Inelastic cross section (mb): " << xsecTot << endl;

	//scale probability distribution to have range (0,1)
	for (int i=0;i<lengths[ind];i++) {
		Xsec[i] /= Xsec[lengths[ind]-1];
	}

	file.close();
	return xsecTot;
}




