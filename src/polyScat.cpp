#include "polyScat.h"
#include <iostream>

polyScat::polyScat(float Zproj0, float thickness)
{
  //thickness is in mg/cm2
  Zproj = Zproj0;
  Zhydrogen = 1.;
  Zcarbon = 6.;
  Nhydrogen =  thickness/1000./16.*2.*6.02e23;
  Ncarbon = thickness/1000./16.*6.02e23;
  Hydrogen = new CMultScat(Zproj,Nhydrogen,Zhydrogen);
  Carbon = new CMultScat(Zproj,Ncarbon,Zcarbon);
}

//Destructor
polyScat::~polyScat(){
  delete Hydrogen;
  delete Carbon;
}

//***********************************************
float polyScat::thetaRMS(float Eproj, float fractionalThick)
{
  float HydrogenRMS = Hydrogen->thetaRMS(Eproj,fractionalThick);
  float CarbonRMS = Carbon->thetaRMS(Eproj,fractionalThick);
  float total = sqrt(pow(HydrogenRMS,2)+pow(CarbonRMS,2));
  return total;
}
