#ifndef _Gobbiarray
#define _Gobbiarray

#include "tele.h"
#include "random.h"
#include <iostream>
using namespace std;


class Gobbiarray
{
 public:
  tele *Tele[4];
  Gobbiarray(float, float, float, float);
  ~Gobbiarray();
  static CRandom ran;


  int hit(float, float, float, float, float, float, float);

  float thetaRecon;
  float phiRecon;
  float xRecon;
  float yRecon;
  int itele;
  float dist;  // Gobbi distance from target center (mm)
  float thick; // Target thickness (mm)
  float length;
  float b;
  float Dactive;
  float thetaCollimator;
  
  int ix;
  int iy;

};


#endif

