#include "gobbi.h"
#include <iostream>
#include <cmath>

using namespace std;

CRandom gobbi::ran;  //random number generator

gobbi::gobbi(double dist0, double b0,bool extra0)
{
  dist = dist0;
  b = b0;//16.;// mm beam axis to frame dimension,
  length = 72.37; // length of Si frame
  Dactive = 64.2; // active dimension of Si
  extra = extra0;

  Tele[0] = new gtele(b+length/2.,length/2.-b,Dactive,length);
  Tele[1] = new gtele(length/2.-b,-b-length/2.,Dactive,length);
  Tele[2] = new gtele(-b - length/2.,-length/2.+b,Dactive,length);
  Tele[3] = new gtele(-length/2.+b,b+length/2.,Dactive,length);


  TeleCsI[0] = new gtele(b+length/2.,length/2.-b,Dactive+5.,length);
  TeleCsI[1] = new gtele(length/2.-b,-b-length/2.,Dactive+5.,length);
  TeleCsI[2] = new gtele(-b - length/2.,-length/2.+b,Dactive+5.,length);
  TeleCsI[3] = new gtele(-length/2.+b,b+length/2.,Dactive+5.,length);

  nTele = 4;

  if (extra)
    {
      Tele[4] = new gtele(3./2.*length-b,-length/2.-b,Dactive,length);
      TeleCsI[4] = new gtele(3./2.*length-b,-length/2.-b,Dactive+5.,length);
      nTele = 5;
    }

  double RadiusCollimator = .75*25.4; //mm
  thetaCollimator = atan(RadiusCollimator/(dist-7.));

}
//*********************************************************
gobbi::~gobbi()
{
  delete Tele[0];
  delete Tele[1];
  delete Tele[2];
  delete Tele[3];
  delete TeleCsI[0];
  delete TeleCsI[1];
  delete TeleCsI[2];
  delete TeleCsI[3];

  if (extra)
    {
      delete Tele[4];
      delete TeleCsI[4];
    }
  
}

//**********************************************************
int gobbi::hit(double theta, double phi, double xtarget, double ytarget, double range, double straggle )
{
  int ihit =0;
  //  if (theta < thetaCollimator) return ihit;
  double x = dist*tan(theta)*cos(phi)+ xtarget;
  double y = dist*tan(theta)*sin(phi)+ ytarget;



  for (int i=0;i<nTele;i++)
    {
      ihit = Tele[i]->hit(x,y);
      if (ihit)
	{
         double r = sqrt(pow(dist,2) + pow(Tele[i]->xRecon,2) 
		     + pow(Tele[i]->yRecon,2));
         thetaHit = acos(dist/r);
         phiHit = atan2(Tele[i]->yRecon, Tele[i]->xRecon);
         itele = i;
         break;
	}

    }
  
  x = (dist+range)*tan(theta)*cos(phi)+ xtarget + ran.Gaus(0.,straggle/sqrt(2.));
  y = (dist+range)*tan(theta)*sin(phi)+ ytarget + ran.Gaus(0.,straggle/sqrt(2.));

  int ihit_CsI = 0;
  for (int i=0;i<nTele;i++)
    {
      ihit_CsI = TeleCsI[i]->hit(x,y);
      if (ihit_CsI)
	{
         itele_CsI = i;
         break;
	}

    }

   if (itele_CsI != itele) return 0;
  if (Tele[itele]->nCsI != TeleCsI[itele]->nCsI) return 0;

  return ihit;
}
//*************************************************
bool gobbi::thruHole(double theta, double phi, double xtarget, double ytarget)
{

  double x = dist*tan(theta)*cos(phi)+ xtarget;
  double y = dist*tan(theta)*sin(phi)+ ytarget;


  for (int i=0;i<4;i++)
    {
      if (Tele[i]->hitAbsorber(x,y)) return false;
    }
  return true;
}
