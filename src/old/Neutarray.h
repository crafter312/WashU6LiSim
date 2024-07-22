#ifndef _Neutarray
#define _Neutarray

#include "random.h"

using namespace std;

class Neutarray
{
 public:
  static CRandom ran;
  double theta_min;
  double theta_max;
  
  Neutarray(double,double);
  ~Neutarray();

  int hit(double,double);
};


#endif

