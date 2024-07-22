#include "Neutarray.h"
#include <cmath>

CRandom Neutarray::ran;
using namespace std;

//created so there can be more complexity added later
Neutarray::Neutarray(double theta_min0, double theta_max0){
  theta_min = theta_min0;
  theta_max = theta_max0;
}
//*********************************************************
Neutarray::~Neutarray(){}

//**********************************************************
int Neutarray::hit(double theta, double energy)
{
  int ihit = 0;
  
  //TODO maybe incorperate xTarget and yTarget, for now keep simple
  
  if (theta > theta_min && theta < theta_max){
    ihit = 1;
  }
  
  //include intrinsic efficiency of p-terphenyl that is ~2cm thick
  if (ihit == 1){
    double prob = 0.27 - 0.27/10.0 * energy;
    double rand = ran.Rndm();
    //if rand>prob takes the event where we don't detect so we switch ihit back to 0
    if (rand > prob){
      ihit = 0;
    }
  }
  
  return ihit;
}

