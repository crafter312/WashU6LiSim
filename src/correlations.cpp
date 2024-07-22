// Addition added to sim.cpp to incorparate differential cross section into the
// selection of angle theta
// Author: Nicolas Dronchi, 07/23/2020
// Modified: Henry Webb, 05/15/2024

#include "random.h"
#include "correlations.h"
#include "constants.h"

using namespace std;

Correlations::Correlations(string* filenamein, string fileelasticin, double E0, double Ex0, double* Ext0s, double* _Xsecs, size_t n)
{
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
	for (int i=1; i<nexits; i++) 
		Xsecs[i] = Xsecs[i-1] + _Xsecs[i];
	for (int i=0; i<nexits; i++)
		Xsecs[i] /= Xsecs[nexits-1];

	// Declare arrays for inelastic input
	lengths   = new int[nexits];
	th_file   = new double*[nexits];
	Xsec_file = new double*[nexits];

	// Some constants and calculations are same for all sampled angles
	setConstants();

	// Read the Fresco fort.201 file containing unscaled elastic cross section data in the center of mass
	readelastic();

	// Read all the supplied Fresco fort.202, fort.203, etc. files containing inelastic cross section data in the center of mass
	for (int i=0; i<nexits; i++)
		readinelastic(filenames[i], i);

	// Randomize theta and phi angles in Center of Mass
	randomAngles();
}

Correlations::~Correlations(){
	delete []th_file;
	delete []Xsec_file;
	delete []th_elastic;
  delete []Xsec_elastic;
}

void Correlations::randomAngles()
{
	phi = 2.*pi*ran.Rndm();

	/**** INELASTIC ****/
	// probability to accept a random theta is dependent on the cross section

	// first pick an exit channel
	float probtr = ran.Rndm();
	int i = 0;
	for (;;)
	{
		if ((Xsecs[i] > probtr) || (i+1 == nexits)) break;
		i++;
	}
	Ext = Exts[i];
	double* Xsec = Xsec_file[i];
	double* th = th_file[i];
	int length = lengths[i];

	// iterate through scaled/integrated cross section until random threshold is met
	probtr = ran.Rndm();
	int j = 0;
	for (;;)
	{
		if ((Xsec[j] > probtr) || (j+1 == length)) break;
		j++;
	}

	// pick random angle between chosen angle and next one
	if (j+1 == length) thetaCM = th[length-1];
	else thetaCM = th[j] + (th[j+1]-th[j])*ran.Rndm();

	/**** ELASTIC ****/

	probtr = ran.Rndm();

	// iterate through scaled/integrated cross section until random threshold is met
	int ii = 0;
	for (;;)
	{
		if ((Xsec_elastic[ii] > probtr) || (ii+1 == lenElastic)) break;
		ii++;
	}
	if (ii+1 == length) thetaElastic = th_elastic[lenElastic-1];
	else thetaElastic = th_elastic[ii] + (th_elastic[ii+1]-th_elastic[ii])*ran.Rndm();
	
	calculateLabAngles();
}

void Correlations::setConstants()
{
	Mp = Mass_7Li/m0;         // mass of projectile Li7
	Mt = Mass_12C/m0;         // mass of target C12
	Mpp = Mass_6Li/m0;        // mass of outgoing projectile Li6
	Mtt = Mass_13C/m0;        // mass of outgoing target C13
	Mred = Mpp*Mtt/(Mpp+Mtt); // reduced mass

	EnergyPA = E/7.0;               // energy per nucleon
	Vbeam = sqrt(2*EnergyPA)*vfact; // beam velocity
	VCM = Vbeam * Mp/(Mp+Mt);       // CM velocity
	VpCM = Vbeam-VCM;               // velocity of projectile in CM frame
	VtCM = VCM;                     // velocity of target in CM frame is CM velocity

	ECMin = 0.5*Mp*pow(VpCM/vfact,2) + 0.5*Mt*pow(VtCM/vfact,2); // kinetic energy of incoming target and projectile in CM frame
	Qrxn = -2.30471; // Q value for 6Li -> alpha + p + n in MeV
}

void Correlations::calculateLabAngles()
{
	// Energy of projectile and target in exit channel
	// Note that the target nucleus excitation energy 
	// changes depending on the exit channel chosen
	double ECMout = ECMin - Qrxn - Exp - Ext;

	double Vrel = sqrt(ECMout*2/Mred)*vfact; // relative velocity between projectile and target in exit channel
	double Vpp = Vrel*Mtt/(Mtt+Mpp);         // center of mass velocity of projectile fragment in exit channel
	double Vtt = Vrel - Vpp;                 // CM velocity of target fragment in exit channel

	// projections of projectile out on x and z in CM frame
	double Vxpp = Vpp*sin(thetaCM);
	double Vzpp = Vpp*cos(thetaCM);

	// projections of projectile out on x and z in Lab frame
	double Vxpplab = Vxpp;
	double Vzpplab = Vzpp + VCM ;

	/**** lab velocity and angle of the projectile ****/
	Vpplab = sqrt((Vxpplab*Vxpplab) + (Vzpplab*Vzpplab));
	thetaLab = acos(Vzpplab/Vpplab);

	// target is always 180 deg away in center of mass frame
	double NangCM = pi - thetaCM;

	// projections of target out on x and z in CM frame
	double Vxtt = Vtt*sin(NangCM);
	double Vztt = Vtt*cos(NangCM);

	// projections of target out on x and z in Lab frame
	double Vxttlab = Vxtt;
	double Vzttlab = Vztt + VCM;

	/**** lab velocity and angle of the target ****/
	Vttlab = sqrt((Vxttlab*Vxttlab) + (Vzttlab*Vzttlab));
	thetaTarg = acos(Vzttlab/Vttlab);
}

void Correlations::readelastic()
{
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

	// Define constants and physical values
	double Mp = Mass_7Li/m0;		 // mass of projectile Li7
	double Mt = Mass_12C/m0;		 // mass of target C12
	double Mred = Mp*Mt/(Mp+Mt); // reduced mass

	double EnergyPA = E/7.0;								 // energy per nucleon
	double Vbeam = sqrt(2*EnergyPA)*vfact; // beam velocity
	double VCM = Vbeam * Mp/(Mp+Mt);			 // velocity of CM
	double VpCM = Vbeam-VCM;							 // velocity of projectile in CM frame
	double VtCM = VCM;										 // velocity of target in CM frame is velocity in CM

	double ECMin = 0.5*Mp*pow(VpCM/vfact,2) + 0.5*Mt*pow(VtCM/vfact,2);
	double ECMout = ECMin;
	
	double Vrel = sqrt(ECMout*2/Mred)*vfact; // relative velocity between projectile and target in exit channel
	double Vp = Vrel*Mt/(Mt+Mp);						 // center of mass velocity of projectile fragment in exit channel

	// Create output file
	fstream ofile;
	ofile.open("transformed_elastic.out", ios::out);
	ofile << "angle,xsec,integral" << endl;

	// Get data from file
	int count = 0;
	double theta_i, theta_rad, Vxplab, Vzplab, Vplab, thetaLab;
	while (getline(file, line))
	{
		if ((line.compare(0,1,"#") == 0) || (line.compare(0,1,"@") == 0) || (line.find("END") != string::npos))
			continue;

		stringstream sstr(line);
		sstr >> theta_i >> Xsec_elastic[count]; // extracts angle and cross section values from a particular input file line

		/**** TRANSFORM ANGLE TO LAB FRAME ****/

		// projections of projectile out on x and z in Lab frame
		theta_rad = theta_i*deg_to_rad;
		Vxplab = Vp*sin(theta_rad);
		Vzplab = Vp*cos(theta_rad) + VCM;

		// lab velocity and angle of the projectile
		Vplab = sqrt((Vxplab*Vxplab) + (Vzplab*Vzplab));
		thetaLab = acos(Vzplab/Vplab);

		th_elastic[count] = thetaLab;
		ofile << to_string(thetaLab) << "," << to_string(Xsec_elastic[count]) << ",";

		// integrate differential cross section for theta > 3 deg.
		if (thetaLab <= 3.*deg_to_rad)
		{
			Xsec_elastic[count] = 0;
			ofile << "0" << endl;
			count++;
			continue;
		}

		Xsec_elastic[count] += Xsec_elastic[count-1];
		ofile << to_string(Xsec_elastic[count]) << endl;
		count++;
	}

	cout << "Elastic cross section (mb): " << Xsec_elastic[lenElastic-2] * ((double) lenElastic) / 180.0 << endl;

	// scale probability distribution to have range (0,1)
	double norm = Xsec_elastic[lenElastic-1];
	for (int i=0; i<lenElastic; i++)
	{
		Xsec_elastic[i] /= norm;
	}

	file.close();
	ofile.close();
}

void Correlations::readinelastic(string filename, int ind)
{
	fstream file;
	cout << "Inelastic Differential Cross Section filename " << filename << endl;
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




