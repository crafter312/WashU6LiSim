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
	CFrame::einstein = einstein; //change einstein for all CFrame objects
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

float CDecay::getErelNewton(CFrame** part)
{

	for (int i=0;i<3;i++) 
	 {
		 plfRecon->v[i] = 0.;
		 for (int j=0;j<Nfrag;j++)
					 plfRecon->v[i] += part[j]->v[i]*part[j]->A;
		 plfRecon->v[i] /= sumA;	
	 }

	plfRecon->getEnergy();

	ErelRecon = 0.;
	 for (int j=0;j<Nfrag;j++)
		 {
			 partCM[j]->velocity = 0.;
			 for (int i=0;i<3;i++)
				 {
					partCM[j]->v[i] = part[j]->v[i] - plfRecon->v[i];
					partCM[j]->velocity += pow(partCM[j]->v[i],2);
				 }
			 partCM[j]->energy = real[j]->A/2.*partCM[j]->velocity/pow(vfact,2);
			 partCM[j]->velocity = sqrt(partCM[j]->velocity);
			 ErelRecon += partCM[j]->energy;
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
	for (int i=0; i<3; i++){
		plfRecon->pc[i] = 0.;
		for (int j=0; j<Nfrag; j++){
			plfRecon->pc[i] += part[j]->pc[i];
		}
	}
	
	//Add up totEnergy of each fragment
	plfRecon->totEnergy = 0;
	for (int j=0;j<Nfrag;j++){
		plfRecon->totEnergy += part[j]->totEnergy;
	}
	
	//start bob changes////////////////////////////////////////////////////
	plfRecon->pcTot = sqrt(pow(plfRecon->pc[0],2) + pow(plfRecon->pc[1],2) + pow(plfRecon->pc[2],2));
	plfRecon->velocity = plfRecon->pcTot/plfRecon->totEnergy*c;
	for (int i=0; i<3; i++){
		plfRecon->v[i] = plfRecon->pc[i]/plfRecon->pcTot*plfRecon->velocity;
	}
	double dv[3];
	for (int i=0; i<3; i++){
		dv[i] = -plfRecon->v[i];
	}
	//end bob changes//////////////////////////////////////////////////////
	
	/* ? had issue with this from Bob changes
	plfRecon->getVelocityFromMom();

	plfRecon->velocity = sqrt(pow(plfRecon->v[0],2)+pow(plfRecon->v[1],2)
				+pow(plfRecon->v[2],2));

	*/
	//cout << "acos(-plfRecon->v[2]/plfRecon->velocity) = acos( " << plfRecon->v[2] << " / " << plfRecon->velocity << " ) " << endl;
	plfRecon->theta = acos(plfRecon->v[2]/plfRecon->velocity);
	double phi = atan2(plfRecon->v[1],plfRecon->v[0]);
	if (phi < 0.) phi += 2.*acos(-1.);
	plfRecon->phi = phi;
	//cout << "theta " << plfRecon->theta*180/acos(-1.0) << endl;

	//transfer real/recon velocities to Center of Mass CFrame
	for (int j=0;j<Nfrag;j++){
		for (int i=0;i<3;i++){
			partCM[j]->v[i] = part[j]->v[i];
		}
	}
	
	ErelRecon = 0.;
	for (int j=0;j<Nfrag;j++){
		//(plfRecon->v) provides reference frame velocity vectors
		partCM[j]->transformVelocity(dv);
		ErelRecon += partCM[j]->getEnergy();
	}

	//emission angle of core - which should be last in the list
	float vv =0.;
	for(int i=0;i<3;i++){
		vv += pow(partCM[Nfrag-1]->v[i],2);
	}
	vv = sqrt(vv);
	cos_thetaH = partCM[Nfrag-1]->v[2]/vv;

	return ErelRecon;
}

float CDecay::getErel_at(CFrame *part1,CFrame *part2){
	//part is either real/recon, both are CFrame objects. plfRecon is also a CFrame
	//add up all the x,y,z momentums of each fragment. The momentum from the
	//decay into fragments will cancle out and give momentum of a particle-like-fragment

	sumA = 0.;
	sumA += part1->A;
	sumA += part2->A;

	plfRecon2 = new CFrame(sumA);

	for (int i=0; i<3; i++){
		plfRecon2->pc[i] = 0.;
		
		plfRecon2->pc[i] += part1->pc[i];
		plfRecon2->pc[i] += part2->pc[i];
		
	}
	
	//Add up totEnergy of each fragment
	plfRecon2->totEnergy = 0;

	plfRecon2->totEnergy += part1->totEnergy;
	plfRecon2->totEnergy += part2->totEnergy;
	

	plfRecon2->pcTot = sqrt(pow(plfRecon2->pc[0],2) + pow(plfRecon2->pc[1],2) + pow(plfRecon2->pc[2],2));
	plfRecon2->velocity = plfRecon2->pcTot/plfRecon2->totEnergy*c;
	for (int i=0; i<3; i++){
		plfRecon2->v[i] = plfRecon2->pc[i]/plfRecon2->pcTot*plfRecon2->velocity;
	}
	double dv[3];
	for (int i=0; i<3; i++){
		dv[i] = -plfRecon2->v[i];
	}

	plfRecon2->theta = acos(plfRecon2->v[2]/plfRecon2->velocity);
	double phi = atan2(plfRecon2->v[1],plfRecon2->v[0]);
	if (phi < 0.) phi += 2.*acos(-1.);
	plfRecon2->phi = phi;


	partCM[0]->A = part1->A;
	partCM[0]->mass = m0*part1->A;
	partCM[1]->A = part2->A;
	partCM[1]->mass = m0*part2->A;

	//transfer real/recon velocities to Center of Mass CFrame
	for (int i=0;i<3;i++){
		partCM[0]->v[i] = part1->v[i];
	}
	for (int i=0;i<3;i++){
		partCM[1]->v[i] = part2->v[i];
	}
	
	ErelRecon = 0.;
	for (int j=0;j<Nfrag;j++){
		//(plfRecon->v) provides reference frame velocity vectors
		partCM[j]->transformVelocity(dv);
		ErelRecon += partCM[j]->getEnergy();
	}

	//emission angle of core - which should be last in the list
	float vv =0.;
	for(int i=0;i<3;i++){
		vv += pow(partCM[Nfrag-1]->v[i],2);
	}
	vv = sqrt(vv);
	cos_thetaH = partCM[Nfrag-1]->v[2]/vv;

	return ErelRecon;
}


//*************************************************************
void CDecay::Mode2Body(double Ex, double gamma, double Q)
{

	//find decay energy, use Breit Wigner distribution is gamma>0
	double ET0 = Ex - Q;
	//use next two lines if simulating background
	//ET;
	//ET = ET0;
	
	if (gamma <= 0.) ET = ET0;
	else 
	{
		for (;;)
		{
			ET = ran.BreitWigner(ET0,gamma);
			//ET = ET0 + gamma*2*(ran.Rndm()-0.5);
			if (ET > 0.0001 && ET < 10.) break;
		}
	} 

	double mu = mass1*mass2/(mass1+mass2);
	double Vrel = sqrt(2.*ET/mu)*vfact;
	double v1	= mass2/(mass1+mass2)*Vrel;
	double v2	= Vrel - v1;


	double gamma1 = 1./sqrt(1.-pow(v1/c,2));
	float pc = gamma1*v1*mass1*m0/c;

	for (;;){
		double E1 = sqrt(pow(pc,2)+pow(mass1*m0,2))- mass1*m0;
		double E2 = sqrt(pow(pc,2)+pow(mass2*m0,2))- mass2*m0;
		double y = E1 + E2	- ET;
		double dE1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2));
		double dE2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2));
		double dy = dE1 + dE2;
		double dpc = - y/dy;
		
		//break if derivative is ~0
		if (fabs(dpc) < .0001) break;
		
		pc += dpc;
	}
	
	v1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2))*c;
	v2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2))*c;


	double theta = acos(2.*ran.Rndm()-1.);
	double phi = 2.*acos(-1.)*ran.Rndm();



	real[0]->v[0] = v1*sin(theta)*cos(phi);
	real[0]->v[1] = v1*sin(theta)*sin(phi);
	real[0]->v[2] = v1*cos(theta);


	real[1]->v[0] = -v2*sin(theta)*cos(phi);
	real[1]->v[1] = -v2*sin(theta)*sin(phi);
	real[1]->v[2] = -v2*cos(theta);

}


//*************************************************************
//
//*************************************************************
void CDecay::ModeLineShapes()
{
	//find decay energy, use Breit Wigner difytribution is gamma>0

	

	for (;;)
	{
		ET = prof->rand_2branches(ran.Rndm(), ran.Rndm());
		if (ET > 0.)
		{
			protonbranch++;
			break;
		}
		else
		{
			neutronbranch++;
		}
	}

	if (ET < 0) { cout << "double check ET selection in ModeLineShapes" << endl; abort();}

	//cout << "ET " << ET << endl;






	double mu = mass1*mass2/(mass1+mass2);
	double Vrel = sqrt(2.*ET/mu)*vfact;
	double v1	= mass2/(mass1+mass2)*Vrel;
	double v2	= Vrel - v1;


	double gamma1 = 1./sqrt(1.-pow(v1/c,2));
	float pc = gamma1*v1*mass1*m0/c;

	for (;;){
		double E1 = sqrt(pow(pc,2)+pow(mass1*m0,2))- mass1*m0;
		double E2 = sqrt(pow(pc,2)+pow(mass2*m0,2))- mass2*m0;
		double y = E1 + E2	- ET;
		double dE1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2));
		double dE2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2));
		double dy = dE1 + dE2;
		double dpc = - y/dy;
		
		//break if derivative is ~0
		if (fabs(dpc) < .0001) break;
		
		pc += dpc;
	}
	
	v1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2))*c;
	v2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2))*c;


	double theta = acos(2.*ran.Rndm()-1.);
	double phi = 2.*acos(-1.)*ran.Rndm();



	real[0]->v[0] = v1*sin(theta)*cos(phi);
	real[0]->v[1] = v1*sin(theta)*sin(phi);
	real[0]->v[2] = v1*cos(theta);


	real[1]->v[0] = -v2*sin(theta)*cos(phi);
	real[1]->v[1] = -v2*sin(theta)*sin(phi);
	real[1]->v[2] = -v2*cos(theta);

}

//**************************************************************
// This function is designed to simulate direct "microcanonical"
// or "phase space" decays in which the momenta of the three
// fragments are chosen randomly from the allowable phase space
//**************************************************************
void CDecay::ModeMicroCanonical(double Ex, double gamma, double Q)
{
	//find decay energy, use Breit Wigner distribution is gamma>0
	double ET0 = Ex - Q;
	//use next two lines if simulating background
	//ET;
	//ET = ET0;
	
	if (gamma <= 0.) ET = ET0;
	else 
	{
		for (;;)
		{
			ET = ran.BreitWigner(ET0,gamma);
			//ET = ET0 + gamma*2*(ran.Rndm()-0.5);
			if (ET > 0.0001 && ET < 10.) break;
		}
	}

	// sample fragment velocities from Gaussian distribution
	valarray <float> vcm(3);
	for (int i=0; i<Nfrag; i++)
	{
		double massn = real[i]->A;
		real[i]->v[0] = ran.Gaus(0.,1.)/massn;
		real[i]->v[1] = ran.Gaus(0.,1.)/massn;
		real[i]->v[2] = ran.Gaus(0.,1.)/massn;
		vcm[0] += real[i]->v[0]*massn;
		vcm[1] += real[i]->v[1]*massn;
	}

	vcm /= sumA;

	// scale velocities to match total energy with desired decay energy
	float testTotal= 0.;
	for (int i=0; i<Nfrag; i++)
	{
		real[i]->velocity = 0.;
		for (int j=0; j<3; j++)
		{
			real[i]->v[j] -= vcm[j];
			real[i]->velocity += real[i]->v[j]*real[i]->v[j];
		}
		real[i]->energy = real[i]->A/2.*real[i]->velocity/(.9784*.9784);
		real[i]->velocity = sqrt(real[i]->velocity);
		testTotal += real[i]->energy;
	}
	float ratio = sqrt(ET/testTotal);
	for (int i=0; i<Nfrag; i++)
	{
		real[i]->velocity *= ratio;
		for (int j=0; j<3; j++) real[i]->v[j] *= ratio;
	}
}

//*************************************************************
void CDecay::Mode2BodyExact(double Ex, double gamma, double Q)
{

	//find decay energy, use Breit Wigner difytribution is gamma>0
	double ET0 = Ex - Q;
	//use next two lines if simulating background
	ET;
	ET = ET0;
 
	double mu = mass1*mass2/(mass1+mass2);
	double Vrel = sqrt(2.*ET/mu)*vfact;
	double v1	= mass2/(mass1+mass2)*Vrel;
	double v2	= Vrel - v1;


	double gamma1 = 1./sqrt(1.-pow(v1/c,2));
	float pc = gamma1*v1*mass1*m0/c;

	for (;;){
		double E1 = sqrt(pow(pc,2)+pow(mass1*m0,2))- mass1*m0;
		double E2 = sqrt(pow(pc,2)+pow(mass2*m0,2))- mass2*m0;
		double y = E1 + E2	- ET;
		double dE1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2));
		double dE2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2));
		double dy = dE1 + dE2;
		double dpc = - y/dy;
		
		//break if derivative is ~0
		if (fabs(dpc) < .0001) break;
		
		pc += dpc;
	}
	
	v1 = pc/sqrt(pow(pc,2)+pow(mass1*m0,2))*c;
	v2 = pc/sqrt(pow(pc,2)+pow(mass2*m0,2))*c;


	double theta = acos(2.*ran.Rndm()-1.);
	double phi = 2.*acos(-1.)*ran.Rndm();



	real[0]->v[0] = v1*sin(theta)*cos(phi);
	real[0]->v[1] = v1*sin(theta)*sin(phi);
	real[0]->v[2] = v1*cos(theta);


	real[1]->v[0] = -v2*sin(theta)*cos(phi);
	real[1]->v[1] = -v2*sin(theta)*sin(phi);
	real[1]->v[2] = -v2*cos(theta);

}
