// Modified by Henry Webb on 24 July 2024

#include "decay.h"

CRandom CDecay::ran;

/**
 * Constructor
\param Nfrag0 gives length of frag array
\param frag0 gives a pointer to an array of Creal objects (6He and p here)
\param einstein0 turns on relativistic calculations when true
 */

CDecay::CDecay(int Nfrag0, CFrag ** frag0, bool einstein0){

	einstein = einstein0;
	Nfrag = Nfrag0;
	frag = frag0;
	
	//transfer info on decay fragments. real/recon are CFrame objects
	for (int i=0;i<Nfrag;i++){
		real[i] = frag[i]->real;
		recon[i] = frag[i]->recon;
	}

	sumA = 0.;
	for (int i=0;i<2;i++) sumA += real[i]->A;
	mass1 = real[0]->A;
	mass2 = real[1]->A;

	plfRecon = new CFrame(sumA);


	for (int i=0;i<Nfrag;i++){
		partCM[i] = new CFrame(real[i]->A);
	}
}


//add on to the constructor becasue some decays will need lineshape input
void CDecay::GenerateProfile(double Ex, double Q)
{
	protonbranch = 0;
	neutronbranch = 0;

	double ET0 = Ex - Q; //10.200-9.973 = 0.227
	double En = Ex - 7.2499; //=2.9501energy between excited state and n+6Li threshold
	double DE = ET0 - En; //=-2.7231


	int z1 = 1; //charge of proton
	int z2 = 2; //charge of 6He
	double mu = 6./7.; //reduced width
	double ac = 4.0; //1.4*(pow(1.,1./3.)+pow(6.,1./3.)); //channel radius
	int l1 = 0; //both are s-wave decays
	int l2 = 0;

	double B_1 = 999.; //boundary conditions
	double B_2 = 999.;

	double WignerLimit = 41.863/mu/pow(ac,2);
	double rwidth2_1 = WignerLimit; //reduced width^2 for proton decay
	double rwidth2_2 = WignerLimit*0.1; //reduced width^2 for neutron decay

	cout << "reduced width^2 proton " << rwidth2_1 << "	reduced width^2 neutron " << rwidth2_2 << endl;

	prof = new profile(ET0, DE, z1, z2, mu, ac, l1, l2, rwidth2_1, rwidth2_2, B_1, B_2);

}







//*********************************************************
	/**
	 * Destructor
	 */
CDecay::~CDecay()
{
	for (int i=0;i<2;i++) delete partCM[i];
	delete plfRecon;
}
//**********************************************************
	/**
	 * returns the reconstructed kinetical energy of the fragmens in their
	 * center-of-mass frame using the real fragment velocities.
	 */
float CDecay::getErelReal()
{
	return getErel(real);
}
//********************************************************
//*********************************************************
	/**
	 * returns the reconstructed kinetical energy of the fragmens in their
	 * center-of-mass frame using the reconstructed or detected
	 *	fragment velocities. 
	 */
float CDecay::getErelRecon()
{
	return getErel(recon);
}
//*********************************************************
float CDecay::getErel(CFrame** part)
{
	if (einstein) return getErelRel(part);
	else return getErelNewton(part);
}
//**********************************************************
	/** 
	 * find the relative kinetic energy of the fragments in their
	 * center-of-mass frame. Non-relativistic version
	 \param part is a pointer to the fragments velocity vectors (real or reconstructed)
	 */

float CDecay::getErelNewton(CFrame** part) {
	plfRecon->SetVelocityComps(0, 0, 0);
	CFrame* tempFrag;
	double A;
	for (int i = 0; i < Nfrag; i++) {
		tempFrag = part[i];
		A = tempFrag->A;
		plfRecon->AddVVec(
			tempFrag->GetVComp(0) * A,
			tempFrag->GetVComp(1) * A,
			tempFrag->GetVComp(2) * A
		);
	}
	plfRecon->ScaleVVec(1. / sumA);
	plfRecon->getEnergy(&einstein);

	ErelRecon = 0.;
	for (int j = 0; j < Nfrag; j++) {
		tempFrag = partCM[j];
		tempFrag->SetVelocityComps(
			part[j]->GetVComp(0) - plfRecon->GetVComp(0),
			part[j]->GetVComp(1) - plfRecon->GetVComp(1),
			part[j]->GetVComp(2) - plfRecon->GetVComp(2)
		);
		tempFrag->CalcVMag();
		tempFrag->SetEnergy(real[j]->A / 2. * tempFrag->GetV2() / vfact2);
		ErelRecon += tempFrag->GetEnergy();
	}

	return ErelRecon;
}

//**********************************************************
	/** 
	 * find the relative kinetic energy of the fragments in their
	 * center-of-mass frame. Relativistic version
	 \param part is a pointer to the fragments velocity vectors (real or reconstructed)
	 */
float CDecay::getErelRel(CFrame **part){
	//part is either real/recon, both are CFrame objects. plfRecon is also a CFrame
	//add up all the x,y,z momentums of each fragment. The momentum from the
	//decay into fragments will cancle out and give momentum of a particle-like-fragment
	plfRecon->SetMomComps(0, 0, 0);
	plfRecon->totEnergy = 0;
	CFrame* tempFrag;
	for (int j = 0; j < Nfrag; j++){
		tempFrag = part[j];
		plfRecon->AddPCVec(tempFrag->GetPCComp(0), tempFrag->GetPCComp(1), tempFrag->GetPCComp(2));
		plfRecon->totEnergy += tempFrag->totEnergy; // add up totEnergy of each fragment
	}
	
	//start bob changes////////////////////////////////////////////////////
	plfRecon->CalcPCMag();
	plfRecon->SetVelocity(plfRecon->GetPC() / plfRecon->totEnergy * c);
	double momVel = plfRecon->GetPC() * plfRecon->GetVelocity();
	plfRecon->SetVelocityComps(
		plfRecon->GetPCComp(0) / momVel,
		plfRecon->GetPCComp(1) / momVel,
		plfRecon->GetPCComp(2) / momVel
	);
	double dv[3] = {
		-plfRecon->GetVComp(0),
		-plfRecon->GetVComp(1),
		-plfRecon->GetVComp(2)
	};
	//end bob changes//////////////////////////////////////////////////////
	
	/* ? had issue with this from Bob changes
	plfRecon->getVelocityFromMom();

	plfRecon->velocity = sqrt(pow(plfRecon->v[0],2)+pow(plfRecon->v[1],2)
				+pow(plfRecon->v[2],2));

	*/
	//cout << "acos(-plfRecon->v[2]/plfRecon->velocity) = acos( " << plfRecon->v[2] << " / " << plfRecon->velocity << " ) " << endl;
	plfRecon->Cart2Sph();
	//cout << "theta " << plfRecon->theta*180/acos(-1.0) << endl;

	//transfer real/recon velocities to Center of Mass CFrame
	//(plfRecon->v) provides reference frame velocity vectors
	ErelRecon = 0.;
	for (int j = 0; j < Nfrag; j++) {
		tempFrag = part[j];
		partCM[j]->SetVelocityComps(tempFrag->GetVComp(0), tempFrag->GetVComp(1), tempFrag->GetVComp(2));
		partCM[j]->transformVelocity(dv, &einstein);
		ErelRecon += partCM[j]->GetEnergy();
	}

	//emission angle of core - which should be last in the list
	//note that the velocity is calculated by the "transfomVelocity" function call
	//TODO: fix ordering of fragments in list
	tempFrag = partCM[Nfrag - 1];
	cos_thetaH = tempFrag->GetVComp(2) / tempFrag->GetVelocity();

	return ErelRecon;
}

float CDecay::getErel_at(CFrame *part1, CFrame *part2) {
	//part is either real/recon, both are CFrame objects. plfRecon is also a CFrame
	//add up all the x,y,z momentums of each fragment. The momentum from the
	//decay into fragments will cancle out and give momentum of a particle-like-fragment

	sumA = part1->A + part2->A;
	plfRecon2 = new CFrame(sumA);
	plfRecon2->SetMomComps(
		part1->GetPCComp(0) + part2->GetPCComp(0),
		part1->GetPCComp(1) + part2->GetPCComp(1),
		part1->GetPCComp(2) + part2->GetPCComp(2)
	);
	
	//Add up totEnergy of each fragment
	plfRecon2->totEnergy = part1->totEnergy + part2->totEnergy;
	plfRecon2->CalcPCMag();
	plfRecon2->SetVelocity(plfRecon2->GetPC() / plfRecon2->totEnergy * c);
	plfRecon2->CalcCartV();
	plfRecon2->ScaleVVec(-1); //TODO: check extra minus sign in CalcCartV function, can potentially remove this line
	
	//TODO: if the minus sign in CalcCartV should stay, then this bit should be removed in addition to the above line
	double dv[3] = {
		-plfRecon2->GetVComp(0),
		-plfRecon2->GetVComp(1),
		-plfRecon2->GetVComp(2)
	};

	plfRecon2->Cart2Sph();

	//transfer real/recon velocities to Center of Mass CFrame
	partCM[0]->A = part1->A;
	partCM[0]->mass = m0 * part1->A;
	partCM[0]->SetVelocityComps(part1->GetVComp(0), part1->GetVComp(1), part1->GetVComp(2));
	partCM[1]->A = part2->A;
	partCM[1]->mass = m0 * part2->A;
	partCM[1]->SetVelocityComps(part2->GetVComp(0), part2->GetVComp(1), part2->GetVComp(2));
	
	//(plfRecon->v) provides reference frame velocity vectors
	//TODO: fix Nfrag reference with regard to the fact that this function only takes two fragments as input
	ErelRecon = 0.;
	for (int j = 0; j < Nfrag; j++){
		partCM[j]->transformVelocity(dv, &einstein);
		ErelRecon += partCM[j]->GetEnergy();
	}

	//emission angle of core - which should be last in the list
	//note that the velocity is calculated by the "transfomVelocity" function call
	cos_thetaH = partCM[Nfrag - 1]->GetVComp(2) / partCM[Nfrag - 1]->GetVelocity();

	return ErelRecon;
}


//*************************************************************
void CDecay::Mode2Body(double Ex, double gamma, double Q) {

	//find decay energy, use Breit Wigner distribution is gamma>0
	double ET0 = Ex - Q;
	//use next two lines if simulating background
	//ET;
	//ET = ET0;
	
	if (gamma <= 0.) ET = ET0;
	else {
		for (;;) {
			ET = ran.BreitWigner(ET0, gamma);
			//ET = ET0 + gamma*2*(ran.Rndm()-0.5);
			if (ET > 0.0001 && ET < 10.) break;
		}
	}

	double mu = mass1 * mass2 / (mass1 + mass2);
	double Vrel = sqrt(2. * ET / mu) * vfact;
	double v1	= mass2 / (mass1 + mass2) * Vrel;
	double v2	= Vrel - v1;
	double gamma1 = 1. / sqrt(1. - (v1 * v1 / c2));
	float pc = gamma1 * v1 * mass1 * m0/ c;

	double denom1, denom2, E1, E2, y, dE1, dE2, dy, dpc;
	for (;;) {
		denom1 = sqrt((pc * pc) + (mass1 * mass1 * m02));
		denom2 = sqrt((pc * pc) + (mass2 * mass2 * m02));
		E1 = denom1 - mass1 * m0;
		E2 = denom2 - mass2 * m0;
		y = E1 + E2 - ET;
		dE1 = pc / denom1;
		dE2 = pc / denom2;
		dy = dE1 + dE2;
		dpc = -y / dy;
		
		//break if derivative is ~0
		if (fabs(dpc) < .0001) break;
		
		pc += dpc;
	}
	
	v1 = pc / sqrt((pc * pc) + (mass1 * mass1 * m02)) * c;
	v2 = pc / sqrt((pc * pc) + (mass2 * mass2 * m02)) * c;
	double theta = acos(2.*ran.Rndm()-1.);
	double phi = 2.*acos(-1.)*ran.Rndm();

	real[0]->SetVelocity(v1);
	real[1]->SetVelocity(-v2);
	real[0]->SetTheta(theta);
	real[1]->SetTheta(theta);
	real[0]->SetPhi(phi);
	real[1]->SetPhi(phi);
	real[0]->Sph2CartV();
	real[1]->Sph2CartV();
}


//*************************************************************
//
//*************************************************************
void CDecay::ModeLineShapes() {
	//find decay energy, use Breit Wigner difytribution is gamma>0

	for (;;) {
		ET = prof->rand_2branches(ran.Rndm(), ran.Rndm());
		if (ET > 0.) {
			protonbranch++;
			break;
		}
		else neutronbranch++;
	}

	if (ET < 0) {
		cout << "double check ET selection in ModeLineShapes" << endl;
		abort();
	}

	//cout << "ET " << ET << endl;

	double mu = mass1 * mass2 / (mass1 + mass2);
	double Vrel = sqrt(2. * ET / mu) * vfact;
	double v1	= mass2 / (mass1 + mass2) * Vrel;
	double v2	= Vrel - v1;
	double gamma1 = 1. / sqrt(1. - (v1 * v1 / c2));
	float pc = gamma1 * v1 * mass1 * m0 / c;

	double denom1, denom2, E1, E2, y, dE1, dE2, dy, dpc;
	for (;;) {
		denom1 = sqrt((pc * pc) + (mass1 * mass1 * m02));
		denom2 = sqrt((pc * pc) + (mass2 * mass2 * m02));
		E1 = denom1 - mass1 * m0;
		E2 = denom2 - mass2 * m0;
		y = E1 + E2 - ET;
		dE1 = pc / denom1;
		dE2 = pc / denom2;
		dy = dE1 + dE2;
		dpc = -y / dy;
		
		//break if derivative is ~0
		if (fabs(dpc) < .0001) break;
		
		pc += dpc;
	}
	
	v1 = pc / sqrt((pc * pc) + (mass1 * mass1 * m02)) * c;
	v2 = pc / sqrt((pc * pc) + (mass2 * mass2 * m02)) * c;
	double theta = acos(2.*ran.Rndm()-1.);
	double phi = 2.*acos(-1.)*ran.Rndm();

	real[0]->SetVelocity(v1);
	real[1]->SetVelocity(-v2);
	real[0]->SetTheta(theta);
	real[1]->SetTheta(theta);
	real[0]->SetPhi(phi);
	real[1]->SetPhi(phi);
	real[0]->Sph2CartV();
	real[1]->Sph2CartV();
}

//**************************************************************
// This function is designed to simulate direct "microcanonical"
// or "phase space" decays in which the momenta of the three
// fragments are chosen randomly from the allowable phase space
//**************************************************************
void CDecay::ModeMicroCanonical(double Ex, double gamma, double Q) {
	double ET0 = Ex - Q;
	//use next two lines if simulating background
	//ET;
	//ET = ET0;

	// Find decay energy, use Breit Wigner distribution if gamma>0
	if (gamma <= 0.) ET = ET0;
	else {
		for (;;) {
			ET = ran.BreitWigner(ET0, gamma);
			//ET = ET0 + gamma*2*(ran.Rndm()-0.5);
			if (ET > 0.0001 && ET < 10.) break;
		}
	}

	// Sample fragment velocities from Gaussian distribution
	double massn;
	valarray<float> vcm(3); // initializes elements to 0
	for (int i = 0; i < Nfrag; i++) {
		massn = real[i]->A;
		real[i]->SetVelocityComps(
			ran.Gaus(0., 1.) / massn,
			ran.Gaus(0., 1.) / massn,
			ran.Gaus(0., 1.) / massn
		);
		vcm[0] += real[i]->GetVComp(0) * massn;
		vcm[1] += real[i]->GetVComp(1) * massn;
	}
	vcm /= sumA;

	// Scale velocities to match total energy with desired decay energy
	CFrame* tempFrag;
	double testTotal = 0.;
	for (int i = 0; i < Nfrag; i++) {
		tempFrag = real[i];
		tempFrag->AddVVec(-vcm[0], -vcm[1], -vcm[2]);
		tempFrag->CalcVMag();
		tempFrag->SetEnergy(tempFrag->A / 2. * tempFrag->GetV2() / vfact2);
		testTotal += tempFrag->GetEnergy();
	}
	double ratio = sqrt(ET / testTotal);
	for (int i = 0; i < Nfrag; i++) {
		tempFrag = real[i];
		tempFrag->SetVelocity(tempFrag->GetVelocity() * ratio);
		tempFrag->ScaleVVec(ratio);
	}
}

//*****************************************
// TODO: Could delete this function?
// Same as case for Mode2Body if gamma <= 0
//*****************************************
void CDecay::Mode2BodyExact(double Ex, double gamma, double Q) {
	double ET0 = Ex - Q;
	//use next two lines if simulating background
	ET;
	ET = ET0;

	double mu = mass1 * mass2 / (mass1 + mass2);
	double Vrel = sqrt(2. * ET / mu) * vfact;
	double v1	= mass2 / (mass1 + mass2) * Vrel;
	double v2	= Vrel - v1;
	double gamma1 = 1. / sqrt(1. - (v1 * v1 / c2));
	float pc = gamma1 * v1 * mass1 * m0 / c;

	double denom1, denom2, E1, E2, y, dE1, dE2, dy, dpc;
	for (;;) {
		denom1 = sqrt((pc * pc) + (mass1 * mass1 * m02));
		denom2 = sqrt((pc * pc) + (mass2 * mass2 * m02));
		E1 = denom1 - mass1 * m0;
		E2 = denom2 - mass2 * m0;
		y = E1 + E2 - ET;
		dE1 = pc / denom1;
		dE2 = pc / denom2;
		dy = dE1 + dE2;
		dpc = -y / dy;
		
		//break if derivative is ~0
		if (fabs(dpc) < .0001) break;
		
		pc += dpc;
	}
	
	v1 = pc / sqrt((pc * pc) + (mass1 * mass1 * m02)) * c;
	v2 = pc / sqrt((pc * pc) + (mass2 * mass2 * m02)) * c;
	double theta = acos(2.*ran.Rndm()-1.);
	double phi = 2.*acos(-1.)*ran.Rndm();

	real[0]->SetVelocity(v1);
	real[1]->SetVelocity(-v2);
	real[0]->SetTheta(theta);
	real[1]->SetTheta(theta);
	real[0]->SetPhi(phi);
	real[1]->SetPhi(phi);
	real[0]->Sph2CartV();
	real[1]->Sph2CartV();
}
