#include "fragneut.h"
#include <cmath>


/**
* Constructor
  \param Z0 is the atomic number of fragment
  \param mass0 is mass of fragment in amu
*/

CFragneut::CFragneut(double mass0){
  
  mass = mass0;
  
  neut = new CFrame(mass);
  
  //pick theta values to describe a ring neutron detector approximated as a ring on a sphere
  theta_min = 65.0/180.0*acos(-1.); //radian
  theta_max = 85.0/180.0*acos(-1.);
  Array = new Neutarray(theta_min, theta_max);

}

//*********************************************************
/**
*Destructor
*/
CFragneut::~CFragneut(){
  delete neut;
  delete Array;
}

//******************************************************************
/**
* Add a velocity vector to the fragments velocity vector.
* Used to transform between reference frames
*/
void CFragneut::setVelocity(double *V){
  neut->v[0] = V[0];
  neut->v[1] = V[1];
  neut->v[2] = V[2];
  neut->getEnergy();
}


//**********************************************************
/**
* logical function determines if a particle is detected
*/
int CFragneut::hit(){
  is_hit = Array->hit(neut->theta, neut->energy);
  return is_hit;
}

