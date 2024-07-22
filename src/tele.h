#ifndef _tele
#define _tele
#include "random.h"
#include <iostream>
using namespace std;

class tele
{
 public:
  float xCenter;
  float yCenter;
  float Dactive;
  float dE_thresh;
  float E_thresh;
  static CRandom ran;


  tele(float,float,float,float,float);
  int hit(float,float,float,float);
  float xRecon;
  float yRecon;
  int ix;
  int iy;
};


#endif
