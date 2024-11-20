// Addition added to sim.cpp to incorparate differential cross section into the
// selection of angle theta
// Author: Nicolas Dronchi, 07/23/2020
// Modified: Henry Webb, 05/15/2024 07/22/2024 11/19/2024

#include "correlations.h"

#include <iostream>

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
	return phi * rad_to_deg;
}

double SampledValues::GetThetaElasticRad() {
	return thetaElastic * rad_to_deg;
}

double SampledValues::GetThetaLabRad() {
	return thetaLab * rad_to_deg;
}

/**********************************************************************************************/

Correlations::Correlations(string* filenamein, string fileelasticin, double E0, double Ex0, double* Ext0s, double* _Xsecs, size_t n, string lossfile_C) {
	CRandom ran;
	filenames = filenamein;
	fileelastic = fileelasticin;

	nexits = n;

	E    = E0;	  // MeV
	Exp  = Ex0;   // MeV
	Exts = Ext0s; // MeV

	// Probability distribution to pick an exit channel using the same logic
	// that's used to create a probability distribution from the differential cross sections
	// NOTE THAT THE LENGTH OF filenamein AND _Xsecs INPUT ARRAYS MUST MATCH nexits!
	Xsecs = new double[nexits];
	Xsecs[0] = _Xsecs[0];
	for (int i = 1; i < nexits; i++) 
		Xsecs[i] = Xsecs[i - 1] + _Xsecs[i];
	for (int i = 0; i < nexits; i++)
		Xsecs[i] /= Xsecs[nexits - 1];

	// Declare arrays for inelastic input
	lengths   = new int[nexits];
	th_file   = new double*[nexits];
	Xsec_file = new double*[nexits];

	// Some constants and calculations are same for all sampled angles
	setConstants();

	// Initialize loss object from file, must be done after setConstants
	ploss_C = new CLoss(lossfile_C, Mp);

	// Read the Fresco fort.201 file containing unscaled elastic cross section data in the center of mass
	readelastic();

	// Read all the supplied Fresco fort.202, fort.203, etc. files containing inelastic cross section data in the center of mass
	for (int i = 0; i < nexits; i++)
		readinelastic(filenames[i], i);
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
	Mp     = Mass_7Li / m0;           // mass of projectile Li7
	Mt     = Mass_12C / m0;           // mass of target C12
	Mpp    = Mass_6Li / m0;           // mass of outgoing projectile Li6
	Mtt    = Mass_13C / m0;           // mass of outgoing target C13
	Mred   = Mpp * Mtt / (Mpp + Mtt); // reduced mass of outgoing projectiles

	// Incoming beam frame
	framep = new CFrame(Mp);
	framep->SetTheta(0.);
	framep->SetPhi(0.);

	// Incoming target frame
	framet = new CFrame(Mt);
	framet->SetTheta(0.);
	framet->SetPhi(0.);
	framet->SetEnergy(0.);
	framet->getVelocityRel();

	// CM frame of incoming projectiles
	frameCMout = new CFrame(Mred);
	framet->SetTheta(0.);
	framet->SetPhi(0.);

	// outgoing parent frame
	framepp = new CFrame(Mpp);       

	Qrxn = mass_7Li + mass_12C - mass_6Li - mass_13C; // Q value for 7Li + 12C -> 6Li + 13C, ~ -2.30478473 MeV
	cout << "n-transfer Q-value: " << Qrxn << endl;
}

void Correlations::calculateLabAngles(double thick) {
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
	//if (sampledValues.Ext == 0.) cout << Vbeam << endl;

	// Energy of projectile and target in exit channel
	// Note that the target nucleus excitation energy 
	// changes depending on the exit channel chosen
	double ECMout = ECMin - Qrxn - Exp - sampledValues.Ext;

	// Set up outgoing CM frame
	frameCMout->SetEnergy(ECMout);
	double VCMout = frameCMout->getVelocityRel();

	double Vrel = sqrt(ECMout * 2 / Mred) * vfact; // relative velocity between projectile and target in exit channel
	double Vpp = Vrel * Mtt / (Mtt + Mpp);         // center of mass velocity of projectile fragment in exit channel
	double Vtt = Vrel - Vpp;                       // CM velocity of target fragment in exit channel

	// projections of projectile out on x and z in CM frame
	double Vxpp = Vpp * sin(thetaCM);
	double Vzpp = Vpp * cos(thetaCM);

	// projections of projectile out on x and z in Lab frame
	double Vxpplab = Vxpp;
	double Vzpplab = Vzpp + VCM ;

	/**** lab velocity and angle of the projectile ****/
	sampledValues.Vpplab = sqrt((Vxpplab * Vxpplab) + (Vzpplab * Vzpplab));
	sampledValues.thetaLab = acos(Vzpplab / sampledValues.Vpplab);
	sampledValues.CalculateCartesian();

	// target is always 180 deg away in center of mass frame
	double NangCM = pi - thetaCM;

	// projections of target out on x and z in CM frame
	double Vxtt = Vtt * sin(NangCM);
	double Vztt = Vtt * cos(NangCM);

	// projections of target out on x and z in Lab frame
	double Vxttlab = Vxtt;
	double Vzttlab = Vztt + VCM;

	/**** lab velocity and angle of the target ****/
	Vttlab = sqrt((Vxttlab * Vxttlab) + (Vzttlab * Vzttlab));
	thetaTarg = acos(Vzttlab / Vttlab);

	// Convert angles to degrees
	sampledValues.phi *= rad_to_deg;
	sampledValues.thetaElastic *= rad_to_deg;
	sampledValues.thetaLab *= rad_to_deg;
}

// Elastic scattering angles are converted to the lab frame when reading the file
void Correlations::readelastic() {
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

	// Calculate constant values
	double MredE = Mp * Mt / (Mp + Mt);             // reduced mass
	double ECMout = ECMin;                          // energy out = energy in (elastic scattering)
	double Vrel = sqrt(ECMout * 2 / MredE) * vfact; // relative velocity between projectile and target in elastic scattering exit channel
	double Vp = Vrel * Mt / (Mt + Mp);						  // center of mass velocity of projectile fragment in exit channel

	// Get data from file
	int count = 0;
	double theta_i, theta_rad, Vrplab, Vzplab, Vplab, thetaLab;
	while (getline(file, line)) {
		if ((line.compare(0,1,"#") == 0) || (line.compare(0,1,"@") == 0) || (line.find("END") != string::npos))
			continue;

		stringstream sstr(line);
		sstr >> theta_i >> Xsec_elastic[count]; // extracts angle and cross section values from a particular input file line

		/**** TRANSFORM ANGLE TO LAB FRAME ****/

		// Projections of projectile out on x and z in Lab frame
		theta_rad = theta_i * deg_to_rad;
		Vrplab    = Vp * sin(theta_rad);
		Vzplab    = Vp * cos(theta_rad) + VCM;

		// Lab velocity and angle of the projectile
		Vplab = sqrt((Vrplab * Vrplab) + (Vzplab * Vzplab));
		thetaLab = acos(Vzplab / Vplab);

		th_elastic[count] = thetaLab;

		// integrate differential cross section for theta > 3 deg.
		Xsec_elastic[count] = (thetaLab > 3. * deg_to_rad) ? Xsec_elastic[count] + Xsec_elastic[count - 1] : 0;
		count++;
	}

	cout << "Elastic cross section (mb): " << Xsec_elastic[lenElastic - 2] * ((double) lenElastic) / 180.0 << endl; // TODO: double check this line

	// scale probability distribution to have range (0, 1)
	double norm = Xsec_elastic[lenElastic - 1];
	for (int i = 0; i < lenElastic; i++)
		Xsec_elastic[i] /= norm;

	file.close();
}

void Correlations::readinelastic(string filename, int ind) {
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
	while (getline(file, line))
	{
		if ((line.compare(0,1,"#") == 0) || (line.compare(0,1,"@") == 0) || (line.find("END") != string::npos))
			continue;

		// store cross section data in arrays
		stringstream sstr(line);
		sstr >> theta_i >> Xsec[count];
		th[count] = theta_i*deg_to_rad;

		// integrate differential cross section
		if (count == 0)
		{
			count++;
			continue;
		}

		Xsec[count] += Xsec[count-1];
		count++;
	}

	//scale probability distribution to have range (0,1)
	for (int i=0;i<lengths[ind];i++)
	{
		Xsec[i] /= Xsec[lengths[ind]-1];
	}

	file.close();
}




