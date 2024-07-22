#include "frame.h"
#include "constants.h"

double const CFrame::m0 = 931.478;
double const CFrame::c=30.;
double const CFrame::pi = acos(-1.);
bool CFrame::einstein = 0;

CFrame::CFrame(double A0)
{
  A = A0;
  mass = m0*A;
  //set initial velocity to 0 (used by Cfragneut but not needed if mode2body is used before)
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;
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
double CFrame::getVelocity()
{
  if (einstein) return getVelocityRel();
  else return getVelocityNewton();
}

//********************************
  /**
   * returns the fragments velocity and calculates its vector (cm/ns)
   * Non-Relativistic version
   */
double CFrame::getVelocityNewton()
{
 velocity =  sqrt(2.*energy/A)*vfact;
 v[2] = velocity*cos(theta);
 v[1] = velocity*sin(theta)*cos(phi);
 v[0] = velocity*sin(theta)*sin(phi);
 return velocity;
}
//*********************************
  /**
   * returns the fragments velocity and calculates it vector (cm/ns)
   * Relativistic version
   */
double CFrame::getVelocityRel()
{

  totEnergy = energy + mass;

  pcTot = sqrt(pow(totEnergy,2)-pow(mass,2));
  pc[2] = pcTot*cos(theta);
  pc[0] = pcTot*sin(theta)*cos(phi);
  pc[1] = pcTot*sin(theta)*sin(phi);

  gamma = totEnergy/mass;
  velocity = c*sqrt(1.-pow(1./gamma,2));
  v[2] = velocity*cos(theta);
  v[0] = velocity*sin(theta)*cos(phi);
  v[1] = velocity*sin(theta)*sin(phi);
  return velocity;
}
//*************************************
double CFrame::getEnergy()
{
  if (einstein) return getEnergyRel();
  else return getEnergyNewton();
}

//***********************************
  /**
   * returns the fragments kinetic energy (MeV)
   * Non-relativistic version
   */
double CFrame::getEnergyNewton()
{
  velocity = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
  energy = 0.5*A*pow(velocity/vfact,2);
  theta = acos(v[2]/velocity);
  phi = atan2(v[1],v[0]);
  if (phi < 0.) phi += 2.*acos(-1.);
  return energy;
}
//***********************************
  /**
   * returns the fragments kinetic energy (MeV)
   * Relativistic version
   */

double CFrame::getEnergyRel()
{
  velocity = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
  gamma = 1./sqrt(1.-pow(velocity/c,2));
  totEnergy = mass*gamma;
  energy = totEnergy - mass;
  pcTot = gamma*velocity*mass/c;
  for (int i=0;i<3;i++) pc[i] = v[i]/velocity*pcTot;
  theta = acos(v[2]/velocity);
  phi = atan2(v[1],v[0]);
  if (phi < 0.) phi += 2.*acos(-1.);
  return energy;

}
//*************************
//*************************************************
void CFrame::getAngle()
{
  velocity = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
  theta = acos(v[2]/velocity);
  phi = atan2(v[1],v[0]);
  if (phi < 0.) phi += 2.*pi;
}
void CFrame::transformVelocity(double *vReference)
{
  if (einstein) transformVelocityRel(vReference);
  else transformVelocityNewton(vReference);
}

//***********************************************
  /**
   * Transforms the velocity vector to a new reference frame
   * Non-Relativistic version
   \param vReference is velocity vector of new reference frame cm/ns
  */
void CFrame::transformVelocityNewton(double * vReference)
{
  for (int i=0;i<3;i++) v[i] += vReference[i];
  getEnergy();
}
//***********************************************
  /**
   * Transforms the velocity vector to a new reference frame
   * Relativistic version
   \param vReference is velocity vector of new reference frame cm/ns
  */
void CFrame::transformVelocityRel(double * vReference)
{
  //find parallel and perpendicular velocities to Vreference
  double vPara[3];
  double vPerp[3];
  double dot = 0.;
  double VVreference = 0.;
  //velocity v is for a specific fragment, vReference is for plf
  //take dot product and magnitude^2 between plf vector and fragment vecotr to
  //be used for a projection
  for (int i=0;i<3;i++){
    dot += v[i]*vReference[i];
    VVreference += pow(vReference[i],2);
  }

  //take projections to find parallel and perpendicular veolocities
  for (int i=0;i<3;i++){
   vPara[i] = dot/VVreference*vReference[i];
   vPerp[i] = v[i] - vPara[i];
  }

  //now transform each
  double vFinalPara[3];
  double vFinalPerp[3];

  double cc = pow(c,2);
  double bb = sqrt(1.-VVreference/cc);

  double vv = 0.;
  for (int i=0;i<3;i++){
    vFinalPara[i] = (vPara[i] + vReference[i])/(1.+dot/cc);
    vFinalPerp[i] = vPerp[i]*bb/(1.+dot/cc);
    v[i] = vFinalPara[i] + vFinalPerp[i];
  }

  getEnergyRel();
}

//***********************************************************
void CFrame::getVelocityFromMom()
{
  if (einstein) getVelocityFromMomRel();
  else getVelocityFromMomNewton();
}

//**************************************
  /**
  * if the momentum vector is known, then gives the velocity vector
  * Relativistic version
  */
void CFrame::getVelocityFromMomRel()
{
  pcTot = 0.;
  for (int i=0;i<3;i++) pcTot += pow(pc[i],2);
  pcTot = sqrt(pcTot);
  totEnergy = sqrt(pow(mass,2) + pow(pcTot,2));
  velocity = pcTot/totEnergy*c;
  for (int i=0;i<3;i++) v[i] = -pc[i]/pcTot*velocity;
  getAngle();
}
//******************************************
void CFrame::getVelocityFromMomNewton()
{
  velocity = 0.;
  for (int i=0;i<3;i++)
    {
     v[i] = -pc[i]/mass*c;
     velocity += pow(v[i],2);
    }
  velocity = sqrt(velocity);
  getAngle();
}
//******************************************
void CFrame::getMomFromVelocity()
{
  velocity = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
  gamma = 1./sqrt(1.-pow(velocity/c,2));
  pcTot = mass*gamma*velocity/c;
  for (int i=0;i<3;i++) pc[i] = v[i]/velocity*pcTot;
  totEnergy = gamma*mass;
}
