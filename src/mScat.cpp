#include "mScat.h"
#include <iostream>

double const CMultScat::a0 = 0.529e-8;
double const CMultScat::pi = acos(-1.0);

/**
* Constructor
\param Zpro0 is particle's atomic number
\param Ztar0 is target's atomic number
\param totalThick0 is targts thickness in mg/cm2
*/

CMultScat::CMultScat(float Zpro0, double totalThick0, float Ztar0){
  Zpro = (double)Zpro0;
  Ztar = (double)Ztar0;
  totalThick = (double)totalThick0;  //atoms per centermeter sqaured

  a = 0.885*a0/sqrt(pow(Zpro,2./3.)+pow(Ztar,2./3.));

  totalTau = pi*pow(a,2)*totalThick;

  factor = 16.26/(Zpro*Ztar)/sqrt(pow(Zpro,2./3.)+pow(Ztar,2./3.))*1000.;
}

/**
* return rms scattering angle 
\param energy is the particles energy in MeV
\param is the fractional thickness of the target through which the particle traverses 
*/

float CMultScat::thetaRMS(float energy, float fractionalThickness)
{
  double tau = (double)fractionalThickness*totalTau;
  double alphaBar = pow(tau,.55);
  double alpha = alphaBar/(double)energy/factor; //half width at half maximum in radians
  
  return (float)(alpha*.851);
}
