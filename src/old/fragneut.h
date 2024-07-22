#ifndef _fragneut
#define _fragneut

#include "Neutarray.h"
#include "frame.h"
#include <string>
#include <iostream>
using std::cout; using std::endl;

class CFragneut
{
 public:
  CFragneut(double);
  ~CFragneut();
  
  double theta_min;
  double theta_max;
  double mass;
  
  Neutarray * Array;
  CFrame *neut;
  
  void setVelocity(double*);
  
  int is_hit;
  int hit();
};

#endif
