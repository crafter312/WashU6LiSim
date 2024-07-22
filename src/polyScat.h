#include "mScat.h"

class polyScat
{
 public:
  polyScat(float,float);
  ~polyScat();
  
  CMultScat* Hydrogen;
  CMultScat* Carbon;
  float Zproj;
  float Zhydrogen;
  float Zcarbon;
  double Nhydrogen;
  double Ncarbon;
  
  float thetaRMS(float,float);
};
