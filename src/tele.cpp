
#include "tele.h"
#include <cmath>


tele::tele( float x0, float y0, float Dactive0, float dE_thresh0, float E_thresh0)
{
  xCenter = x0;
  yCenter = y0;
  Dactive = Dactive0;
  dE_thresh = dE_thresh0;
  E_thresh = E_thresh0;

}

//************************************************
int tele::hit(float x, float y, float dE, float E)
{
  int ihit =0;

  if (fabs(x - xCenter) > Dactive/2.0) return ihit;
  if (fabs(y - yCenter) > Dactive/2.0) return ihit;

  //check energy thresholds to determine if it hits
  if (dE < dE_thresh || E < E_thresh) return ihit;


  ihit = 1;

  ix = (int)((x - xCenter - Dactive/2.)/Dactive*32.);
  iy = (int)((y - yCenter - Dactive/2.)/Dactive*32.);

  xRecon = ((float)ix+ran.Rndm())/32.*Dactive + Dactive/2. + xCenter;
  yRecon = ((float)iy+ran.Rndm())/32.*Dactive + Dactive/2. + yCenter;

  return ihit;
}
