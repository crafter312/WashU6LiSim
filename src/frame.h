#ifndef frame_
#define frame_
#include <cmath>


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
  static double const m0; //!< nucleon mass in MeV
  static double const c; //!< speed of light (cm/ns)
  static double const pi; //!< 3.14159....

  double getVelocity();
  double getVelocityNewton();
  double getVelocityRel();
  double getEnergy();
  double getEnergyFromMomentum();
  double getEnergyNewton();
  double getEnergyRel();
  void getAngle();
  void transformVelocity(double*);
  void transformVelocityNewton(double*);
  void transformVelocityRel(double*);
  void getVelocityFromMom();
  void getVelocityFromMomNewton();
  void getVelocityFromMomRel();
  void getMomFromVelocity();
  CFrame(double);
  static bool  einstein;
};


#endif
