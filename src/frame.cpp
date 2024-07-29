// Modified by Henry Webb on 24 July 2024

#include "frame.h"

KinematicValues::KinematicValues() {
	v[0] = 0; // Set initial velocity to 0 (used by Cfragneut but not needed if mode2body is used before)
  v[1] = 0;
  v[2] = 0;
	pc[0] = NAN;
  pc[1] = NAN;
  pc[2] = NAN;
}

KinematicValues::~KinematicValues() {}

void KinematicValues::Clear() {
	energy = NAN;
	velocity = NAN;
	pcTot = NAN;
	v[0] = 0; // Set initial velocity to 0 (used by Cfragneut but not needed if mode2body is used before)
  v[1] = 0;
  v[2] = 0;
	pc[0] = NAN;
  pc[1] = NAN;
  pc[2] = NAN;
	theta = NAN;
	phi = NAN;
	x = NAN;
	y = NAN;
}

// Converts V from spherical to cartesian coordinates
void KinematicValues::Sph2CartV() {
	v[0] = velocity * sin(theta) * cos(phi);
	v[1] = velocity * sin(theta) * sin(phi);
	v[2] = velocity * cos(theta);
}

// Converts PC from spherical to cartesian coordinates
void KinematicValues::Sph2CartPC() {
  pc[0] = pcTot*sin(theta)*cos(phi);
  pc[1] = pcTot*sin(theta)*sin(phi);
	pc[2] = pcTot*cos(theta);
}

// Calculates spherical angles from V components
void KinematicValues::Cart2Sph() {
	theta = acos(v[2] / velocity);
	phi = atan2(v[1], v[0]);
	phi += (phi < 0) * 2. * pi;
}

// Calculates total V from components
void KinematicValues::CalcVMag() {
	velocity = sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
}

// Calculates total PC from components
void KinematicValues::CalcPCMag() {
	pcTot = sqrt((pc[0] * pc[0]) + (pc[1] * pc[1]) + (pc[2] * pc[2]));
}

// Calculates V components from total PC, total V, and PC components
void KinematicValues::CalcCartV() {
	v[0] = pc[0] / pcTot * velocity;
	v[1] = pc[1] / pcTot * velocity;
	v[2] = pc[2] / pcTot * velocity;
}

// Calculates PC components from total V, total PC, and V components 
void KinematicValues::CalcCartPC() {
	pc[0] = v[0] / velocity * pcTot;
	pc[1] = v[1] / velocity * pcTot;
	pc[2] = v[2] / velocity * pcTot;
}

// Adds supplied vector to V vector
void KinematicValues::AddVVec(double* vec) {
	v[0] += vec[0];
	v[1] += vec[1];
	v[2] += vec[2];
}

// Adds supplied vector to PC vector
void KinematicValues::AddPCVec(double* vec) {
	pc[0] += vec[0];
	pc[1] += vec[1];
	pc[2] += vec[2];
}

// Returns V squared
double KinematicValues::GetV2() {
	return velocity * velocity;
}

// Returns PC squared
double KinematicValues::GetPC2() {
	return pcTot * pcTot;
}

// Returns dot product of V vector with supplied vector
double KinematicValues::GetVDot(double* vec) {
	return (v[0] * vec[0]) + (v[1] * vec[1]) + (v[2] * vec[2]);
}

/********************************************************/

CFrame::CFrame(double A0) {
  A = A0;
  mass = m0 * A;
	mass2 = mass * mass;
}
//******************************************************
/*
double CFrame::getEnergyFromMomentum()
{
  pcTot = sqrt(pow(pc[0],2) + pow(pc[1],2) + pow(pc[2],2));
  theta = acos(pc[2]/pcTot);
  phi = atan2(pc[1],pc[0]);
  totEnergy = sqrt(pow(mass,2)+pow(pcTot,2));
  gamma = totEnergy/mass;
  velocity = c*sqrt(1.-pow(1./gamma,2));
  v[2] = velocity*cos(theta);
  v[0] = velocity*sin(theta)*cos(phi);
  v[1] = velocity*sin(theta)*sin(phi);
  energy = totEnergy - mass;
  return energy;
}
*/

//******************************************************
double CFrame::getVelocity(bool* einstein) {
  if (*einstein) return getVelocityRel();
  else return getVelocityNewton();
}

//********************************
  /**
   * returns the fragments velocity and calculates its vector (cm/ns)
   * Non-Relativistic version
   */
double CFrame::getVelocityNewton() {
	kinematicValues.velocity = sqrt(2. * kinematicValues.energy / A) * vfact;
	kinematicValues.Sph2CartV();
	return kinematicValues.velocity;
}
//*********************************
  /**
   * returns the fragments velocity and calculates it vector (cm/ns)
   * Relativistic version
   */
double CFrame::getVelocityRel() {
  totEnergy = kinematicValues.energy + mass;

  kinematicValues.pcTot = sqrt((totEnergy * totEnergy) - mass2);
  kinematicValues.Sph2CartPC();

  gamma = totEnergy / mass;
  kinematicValues.velocity = c * sqrt(1. - (1. / (gamma * gamma)));
  kinematicValues.Sph2CartV();
  return kinematicValues.velocity;
}
//*************************************
double CFrame::getEnergy(bool* einstein) {
  if (*einstein) return getEnergyRel();
  else return getEnergyNewton();
}

//***********************************
  /**
   * returns the fragments kinetic energy (MeV)
   * Non-relativistic version
   */
double CFrame::getEnergyNewton() {
  kinematicValues.CalcVMag();
  kinematicValues.energy = 0.5 * A * kinematicValues.GetV2() / vfact2;
  kinematicValues.Cart2Sph();
  return kinematicValues.energy;
}
//***********************************
  /**
   * returns the fragments kinetic energy (MeV)
   * Relativistic version
   */

double CFrame::getEnergyRel() {
  kinematicValues.CalcVMag();
  gamma = 1. / sqrt(1. - (kinematicValues.GetV2() / c2));
  totEnergy = mass * gamma;
  kinematicValues.energy = totEnergy - mass;
  kinematicValues.pcTot = gamma * kinematicValues.velocity * mass / c;
  kinematicValues.CalcCartPC();
  kinematicValues.Cart2Sph();
  return kinematicValues.energy;
}
//*************************
//*************************************************
void CFrame::transformVelocity(double* vReference, bool* einstein) {
  if (*einstein) transformVelocityRel(vReference);
  else transformVelocityNewton(vReference);
}

//***********************************************
  /**
   * Transforms the velocity vector to a new reference frame
   * Non-Relativistic version
   \param vReference is velocity vector of new reference frame cm/ns
  */
void CFrame::transformVelocityNewton(double* vReference) {
  kinematicValues.AddVVec(vReference);
  getEnergyNewton();
}
//***********************************************
  /**
   * Transforms the velocity vector to a new reference frame
   * Relativistic version
   \param vReference is velocity vector of new reference frame cm/ns
  */
void CFrame::transformVelocityRel(double* vReference) {
  //find parallel and perpendicular velocities to Vreference
  double vPara[3];
  double vPerp[3];

	//velocity v is for a specific fragment, vReference is for plf
  //take dot product and magnitude^2 between plf vector and fragment vecotr to
  //be used for a projection
  double dot = kinematicValues.GetVDot(vReference);
  double VVreference = (vReference[0] * vReference[0]) + (vReference[1] * vReference[1]) + (vReference[2] * vReference[2]);

  //take projections to find parallel and perpendicular veolocities
  for (int i = 0; i < 3; i++) {
   vPara[i] = dot / VVreference * vReference[i];
   vPerp[i] = kinematicValues.v[i] - vPara[i];
  }

  //now transform each
  double vFinalPara[3];
  double vFinalPerp[3];

  double bb = sqrt(1. - VVreference / c2);

  double vv = 0.;
  for (int i = 0; i < 3; i++) {
		vFinalPara[i] = (vPara[i] + vReference[i]) / (1. + dot / c2);
		vFinalPerp[i] = vPerp[i] * bb / (1. + dot / c2);
		kinematicValues.v[i] = vFinalPara[i] + vFinalPerp[i];
	}

  getEnergyRel();
}

//***********************************************************
void CFrame::getVelocityFromMom(bool* einstein) {
  if (*einstein) getVelocityFromMomRel();
  else getVelocityFromMomNewton();
}

//**************************************
  /**
  * if the momentum vector is known, then gives the velocity vector
  * Relativistic version
  */
void CFrame::getVelocityFromMomRel() {
  kinematicValues.CalcPCMag();
  totEnergy = sqrt(mass2 + kinematicValues.GetPC2());
  kinematicValues.velocity = kinematicValues.pcTot / totEnergy * c;
	kinematicValues.CalcCartV();
  kinematicValues.CalcVMag();
  kinematicValues.Cart2Sph();
}
//******************************************
void CFrame::getVelocityFromMomNewton() {
	kinematicValues.v[0] = -kinematicValues.pc[0] / mass * c; //TODO: double check extra minus sign here
	kinematicValues.v[1] = -kinematicValues.pc[1] / mass * c;
	kinematicValues.v[2] = -kinematicValues.pc[2] / mass * c;
  kinematicValues.CalcVMag();
  kinematicValues.Cart2Sph();
}
//******************************************
void CFrame::getMomFromVelocity() {
  kinematicValues.CalcVMag();
  gamma = 1. / sqrt(1. - (kinematicValues.GetV2() / c2));
  kinematicValues.pcTot = mass * gamma * kinematicValues.velocity / c;
	kinematicValues.CalcCartPC();
  totEnergy = gamma * mass;
}

void CFrame::AddVVec(double v1x, double v1y, double v1z) {
	double v1[3] = { v1x, v1y, v1z };
	kinematicValues.AddVVec(v1);
}

void CFrame::AddPCVec(double pc1x, double pc1y, double pc1z) {
	double pc1[3] = { pc1x, pc1y, pc1z };
	kinematicValues.AddPCVec(pc1);
}

void CFrame::ScaleVVec(double s) {
	kinematicValues.v[0] *= s;
	kinematicValues.v[1] *= s;
	kinematicValues.v[2] *= s;
}

void CFrame::CalcVMag() {
	kinematicValues.CalcVMag();
}

void CFrame::CalcPCMag() {
	kinematicValues.CalcPCMag();
}

void CFrame::Cart2Sph() {
	kinematicValues.Cart2Sph();
}

void CFrame::CalcCartV() {
	kinematicValues.CalcCartV();
}

void CFrame::Sph2CartV() {
	kinematicValues.Sph2CartV();
}

double CFrame::GetV2() {
	return kinematicValues.GetV2();
}

void CFrame::SetEnergy(double E) {
	kinematicValues.energy = E;
}

void CFrame::SetVelocity(double v) {
	kinematicValues.velocity = v;
}

void CFrame::SetVelocityComps(double vx, double vy, double vz) {
	kinematicValues.v[0] = vx;
	kinematicValues.v[1] = vy;
	kinematicValues.v[2] = vz;
}

void CFrame::SetMomComps(double pcx, double pcy, double pcz) {
	kinematicValues.pc[0] = pcx;
	kinematicValues.pc[1] = pcy;
	kinematicValues.pc[2] = pcz;
}

void CFrame::SetTheta(double th) {
	kinematicValues.theta = th;
}

void CFrame::SetPhi(double ph) {
	kinematicValues.phi = ph;
}

void CFrame::SetXY(double x, double y) {
	kinematicValues.x = x;
	kinematicValues.y = y;
}

double CFrame::GetEnergy() {
	return kinematicValues.energy;
}

double CFrame::GetVelocity() {
	return kinematicValues.velocity;
}

double CFrame::GetPC() {
	return kinematicValues.pcTot;
}

double CFrame::GetVComp(int i) {
	return kinematicValues.v[i];
}

double CFrame::GetPCComp(int i) {
	return kinematicValues.pc[i];
}

double CFrame::GetTheta() {
	return kinematicValues.theta;
}

double CFrame::GetPhi() {
	return kinematicValues.phi;
}

double CFrame::GetX() {
	return kinematicValues.x;
}

double CFrame::GetY() {
	return kinematicValues.y;
}
