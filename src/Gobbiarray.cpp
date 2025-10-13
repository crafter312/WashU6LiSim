#include "constants.h"
#include "Gobbiarray.h"
#include <cmath>

using namespace std;

// dist0 - Gobbi distance from target center (mm)
// th - target thickness (mm)
// b0 - Gobbi frame inner dimension (mm)
// radColl - Gobbi collimator outer radius (mm)
Gobbiarray::Gobbiarray(float dist0, float th, float b0, float radColl)
{
  dist = dist0;
  thick = th;
  b = b0;// mm beam axis to frame dimension,
  length = 72.37; // length of Si frame
  Dactive = 64.2; // active dimension of Si

                    //xCenter, yCenter, Dactive,       dE-threshold, E-threshold
  Tele[0] = new tele(b+length/2.,length/2.-b,Dactive,     0.2, 0.5);
  Tele[1] = new tele(length/2.-b,-b-length/2.,Dactive,    0.2, 0.5);
  Tele[2] = new tele(-b - length/2.,-length/2.+b,Dactive, 0.2, 0.5);
  Tele[3] = new tele(-length/2.+b,b+length/2.,Dactive,    0.2, 0.5);

  thetaCollimator = atan(radColl/(dist-7.));
	cout << "Collimator outer angle (degrees): " << thetaCollimator*rad_to_deg << endl;
}
//*********************************************************
Gobbiarray::~Gobbiarray()
{
  delete Tele[0];
  delete Tele[1];
  delete Tele[2];
  delete Tele[3];
}
//**********************************************************
// inthick - the percentage of the target thickness from the upstream side at which the reaction happened
int Gobbiarray::hit(float theta, float phi, float inthick, float xTarget, float yTarget, float dE, float E)
{
  int ihit = 0;
  
  //Hit Beam blocker donut
  if (theta < thetaCollimator) { return ihit;}
  
  float distTot = dist - (0.5 * inthick * thick);
  float x = distTot*tan(theta)*cos(phi) + xTarget;
  float y = distTot*tan(theta)*sin(phi) + yTarget;

  
  //search each quadrant of gobbi array
  for (int i=0;i<4;i++){
    ihit = Tele[i]->hit(x,y,dE,E);
    float r = sqrt(pow(dist,2) + pow(Tele[i]->xRecon,2) 
       + pow(Tele[i]->yRecon,2));

    xRecon = Tele[i]->xRecon;
    yRecon = Tele[i]->yRecon;
    thetaRecon = acos(dist/r);


    phiRecon = atan2(Tele[i]->yRecon, Tele[i]->xRecon);
    itele = i;
    //save the strip number that was hit
    ix = Tele[i]->ix;
    iy = Tele[i]->iy;
    
    if (ihit) break;
  }
  
  return ihit;
}

