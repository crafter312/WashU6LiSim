// Modified by Henry Webb on 24 July 2024

#ifndef frame_
#define frame_

#include <cmath>

#include "constants.h"


/**
 *!\brief Relativistic and Non-Relativistic stuff
 *
 * give velocity, energies, transforms to new reference frame
 * using either non-relativistic or Relativistics equations
 */

class CFrame
{
 public:
  double energy;  //!< fragments energy in MeV
  double velocity; //!< fragments velocity in cm/ns
  double pcTot;  //!< total momentum*c of fragment
  double v[3];  //!< velocity vector of fragment
  double pc[3]; //!< momentum*c vector of fragment MeV
  double mass;  //!< rest mass of fragment in MeV
  double A;
  double gamma; //!< relativistic gamma factor
  double theta; //!< theta agle of fragment in radians
  double phi;  //!< phi angle
  double x;
  double y;
  double totEnergy; //!< rest mass plus kinetic energy

  double getVelocity(bool*);
  double getVelocityNewton();
  double getVelocityRel();
  double getEnergy(bool*);
  //double getEnergyFromMomentum();
  double getEnergyNewton();
  double getEnergyRel();
  void getAngle();
  void transformVelocity(double*, bool*);
  void transformVelocityNewton(double*, bool*);
  void transformVelocityRel(double*);
  void getVelocityFromMom(bool*);
  void getVelocityFromMomNewton();
  void getVelocityFromMomRel();
  void getMomFromVelocity();
  CFrame(double);
};


#endif
