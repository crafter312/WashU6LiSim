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
  Gobbiarray(float);
  ~Gobbiarray();
  static CRandom ran;


  int hit(float,float,float,float,float,float);

  float thetaRecon;
  float phiRecon;
  float xRecon;
  float yRecon;
  int itele;
  float dist;
  float length;
  float b;
  float Dactive;
  float thetaCollimator;
  
  int ix;
  int iy;

};


#endif

