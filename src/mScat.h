#ifndef scat_
#define scat_

#include <cmath>
#include <iostream>

using namespace std;

/**
 *!\brief Multiple scattering of a particle in the target
 */

class CMultScat
{
 public:
 
  double a; //!<screening parameter in fermi's
  double totalThick; //!<target thickness in mg/cm2
  double totalTau;
  double factor;
  static double const a0;
  static double const pi; //!< 3.14159......

  double Zpro;  //!< atomic number of projectile
  double Ztar;  //!< atomic number of target

  CMultScat(float,double,float); 
  float thetaRMS(float,float);
};

#endif
