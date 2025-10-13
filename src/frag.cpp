#include "frag.h"
#include <cmath>


float const CFrag::pi=acos(-1.);
CRandom CFrag::ran;

/**
 * Constructor
 \param Z0 is the atomic number of fragment
 \param mass0 is mass of fragment in amu
 \param filename specifies the name of the file with energy loss information
 \param CsI_res0 fractional energy resolution of this fragment in CsI
 \param thickness is target thickness in mg/cm2 
 \param scaleSmallAngle - scale the predicted width for small angle scattering
 \param	
*/

CFrag::CFrag(float Z0, float mass0, string lossfile_C, string lossfile_Si, float CsI_res0,
	float thickness, shared_ptr<Gobbiarray> array0, float scaleSmallAngle0, bool einstein0, bool useRealP0)
	: Array(array0) {			 
	cout << "Fragment of Z = " << Z0 << ", mass = " << mass0 << endl;	
	Z = Z0;
	mass = mass0;
	//ifstream ifile;
	//ifile.open("~/li6sim/src/teles.dat");

	extra = false;
	loss_C = new CLoss(lossfile_C,mass);
	loss_Si = new CLoss(lossfile_Si,mass);
	CsI_res = CsI_res0;

	//molar mass of CD2 target is 18g/mol
	//float thick = thickness/1000./18.*6.02e23; // atoms/cm2
	
	pScat = new polyScat(Z, thickness);

	real = new CFrame(mass);
	recon = new CFrame(mass);

	scaleSmallAngle = scaleSmallAngle0;
	einstein = einstein0;
	useRealP = useRealP0;
	stopped = false;
}


//*********************************************************
/**
*Destructor
*/
CFrag::~CFrag() {
	if (recon != real) delete recon;
	delete real;
}

//**********************************************************
/**
* logical function determines if a particle is detected
* inthick - the percentage of the target thickness from the upstream side at which the reaction happened
*/
int CFrag::hit(float inthick, float xTarget, float yTarget) {
	//stopped is calculated in targetinteraction, determines if fragments is stopped in the target or in the DeltaE
	if (stopped) return -1;

	is_hit = Array->hit(real->GetTheta(), real->GetPhi(), inthick, xTarget, yTarget, DeltaEnergy, FrontEnergy) ;

	if (is_hit) {
		recon->SetTheta(Array->thetaRecon);
		recon->SetPhi(Array->phiRecon);
		recon->SetXY(Array->xRecon, Array->yRecon);
		recon->SetEnergy(real->GetEnergy() * (1 + CsI_res * ran.Gaus(0.,1.)));
		recon->getVelocity(&einstein);
	}

	return is_hit;
}

//**********************************************************
	/**
	 * return what strips were hit in the int arrays passed in
	 */
void CFrag::getStripHit(int *stripx, int *stripy, int index) {
	stripx[index] = Array->ix;
	stripy[index] = Array->iy;
}


//******************************************************************
/**
* Add a velocity vector to the fragments velocity vector.
* Used to transform between reference frames
*/
void CFrag::AddVelocity(double *Vplf){ 
	real->transformVelocity(Vplf, &einstein);
}
//****************************************************************
/**
* returns the energy after the fragment has exited the target
\param thick is the thickness of target that the particle has to traverse (mg/cm2)
*/
float CFrag::Eloss(float thick){
	if (real->GetEnergy() <= 0.) return 0.;
	real->SetEnergy(loss_C->getEout(real->GetEnergy(), thick));
	return real->GetEnergy();
}
//*******************************************************************
/**
* corrects energy of a detected particle for the energy loss
* in the target.
\param thick is the thickness of target material though which the particle passed (mg/cm2)
*/
float CFrag::Egain(float thick) {
	if (useRealP) return EgainHelper(thick, real);
	return EgainHelper(thick, recon);
}

float CFrag::EgainHelper(float thick, CFrame* frame) {
	if (thick > 0.)
		frame->SetEnergy(loss_C->getEin(frame->GetEnergy(), thick / cos(frame->GetTheta())));

	frame->getVelocity(&einstein);
	return frame->GetEnergy();
} 
//***********************************************
//include multiple scattering
/**
* Monte Carlo choice of small angle scattering due to passage through the target
\param fractionalThick is the fractional thick of the target through which the particle passed
*/
void CFrag::MultiScat(float fractionalThick){
//float Zproj, float Eproj, float fractionalThick
	float thetaRMS = pScat->thetaRMS(real->GetEnergy(), fractionalThick);
	//cout << "thetaRMS " << thetaRMS*180/acos(-1.0) << endl;
	float sigma = thetaRMS/sqrt(2.)*scaleSmallAngle;
	//cout << "thetaRMS= " << thetaRMS << endl;
	float deltaTheta = sqrt(2.)*sigma*sqrt(-log(ran.Rndm()));
	//cout << "deltaTheta= " << deltaTheta << endl;
	float deltaPhi = 2.*pi*ran.Rndm();
	//cout << "delta Phi= " << deltaPhi << endl;

	float x = sin(deltaTheta)*cos(deltaPhi);
	float y = sin(deltaTheta)*sin(deltaPhi);
	float z = cos(deltaTheta);



	//rotate in z-x plane by theta
	double theta = real->GetTheta();
	float xx = x*cos(theta) + z*sin(theta);
	float yy = y;
	float zz = z*cos(theta) - x*sin(theta);


	//rotate in x-y plane
	double phi = real->GetPhi();
	float xxx = xx*cos(phi) - yy*sin(phi);
	float yyy = yy*cos(phi) + xx*sin(phi);
	float zzz = zz;


	float thetaNew = acos(zzz);
	float phiNew = atan2(yyy,xxx);


	real->SetTheta(thetaNew);
	real->SetPhi(phiNew);
}
//*********************
/**
* accounts for multiscattering and energy loss in the target
\param dthick is thickness of target though the particle must pass (mg/cm2)
\param thickness is total target thickness (mg/cm2)
	 */
bool CFrag::targetInteraction(float dthick, float thickness) {

	stopped = 0;
	//protect against random zero dthick
	if (dthick == 0.) return stopped;
	//after scattering, takes path through target at angle
	float thick = dthick / cos(real->GetTheta());
	//if (ran.Rndm() < 0.25) thick = thick*1.5;
	
	Eloss(thick); //returns the energy after the fragment has exited the target
	
	//check if fragment got stuck within the target
	if (real->GetEnergy() <= 0.) { 
		stopped = 1;
		return stopped;
	}

	//changes theta and phi based off of RMStheta of scattering
	MultiScat(thick / thickness);

	return stopped;
}

//*********************
/**
* gets the energy deposited in the DeltaE and E detector needed for checking thresholds
*/
bool CFrag::SiliconInteraction() {
	//no need to look at particles stopped in the target
	if (stopped) return stopped;

	//after scattering, takes path through Si at angle
	float dE_thick = 14.85 / cos(real->GetTheta()); //14.85 mg/cm^2 is 0.064mm thick Si
	
	FrontEnergy = loss_Si->getEout(real->GetEnergy(), dE_thick);
	DeltaEnergy = real->GetEnergy() - FrontEnergy;

	//check if particle was stopped in the dE detector, won't have PID
	if (FrontEnergy <= 0.) {
		DeltaEnergy = -1;
		stopped = 1;
		return stopped;
	}

	return stopped;
}


