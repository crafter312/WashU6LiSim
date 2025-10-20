#include "loss.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

/**
 * constructor
\param filename is name of file containing energy loss tables of a particulat particle
\param mass0 is mass of particle in amu
*/

CLoss::CLoss(string filename, float mass0){
  mass = mass0;
  
  string longfilename = string(LOSSPATH) + "/" + filename;
  ifstream File(longfilename.c_str());
  if (File.is_open() != 1){
    cout << " could not open loss file " << longfilename;
    return;
  }
  
  //get descriptive header of the file
  char line[120];
  File.getline(line,120);
  cout << line << endl;

  File >> N;
  Ein = new float [N];
  dedx = new float [N];
  
  //input data from file into arrays
  for (int i=0;i<N;i++){
    File >> Ein[i] >> dedx[i];
    //cout << Ein[i] << " " << dedx[i] << endl;
  }
}

//****************************************************************
/**
* destructor
*/
CLoss::~CLoss(){
  delete [] Ein;
  delete [] dedx;
}

//*****************************************************************
/*
* returns the value of DeDx interpolated from table
\param energy is energy of particle in MeV
*/
float CLoss::getDedx(float energy){
  //linear interpolation
  int istart = 0;
  float epa = energy/mass;
  if (epa < Ein[0]) istart = 0; //check if before first energy in table
  else if (epa > Ein[N-1]) istart =  N-1; //check if after last energy in table
  else{ //find the index in table that has energy just less than required
    istart = 1;
    for (;;){
      if (epa < Ein[istart]) break;
      istart++;
	  }
    istart--;
  }
  //interpolate between energy just less than and just more than
  float de = (epa-Ein[istart])/(Ein[istart+1]-Ein[istart])
    *(dedx[istart+1]-dedx[istart]) + dedx[istart];
  
  return de;
}

//********************************************************************
/**
* returns the residual energy of particle after passage through absorber
\param energy is initial energy of particle in MeV
\param thick is the thickness of the absorber in mg/cm2
*/
float CLoss::getEout(float energy, float thick)
{
  float dthick = 0.01;
  float de;
  float Eout= energy;
  
  //take small steps of thickness to decrease energy of fragment
  while(thick > 0.0){
    float thickness = min(thick,dthick);
    de = getDedx(Eout); //get interp value
    Eout -= de*thickness;
    if (Eout < 0. || Eout != Eout) return -1.;
    thick -= dthick;
  }
  
  /* old loop
  for(;;){
     float thickness = min(thick,dthick);
     de = getDedx(Eout);
     Eout -= de*thickness;
     if (Eout <0. || Eout != Eout) return -1.;
     if (thickness == thick) break;
     thick -= dthick;
   }
   */
   return Eout;
}
//********************************************************************
  /**
   * Determined the inital energy of a particle before entering an 
   * absorber given its residueal
\param energy is the residual energy of the particle
\param thick is the thickness of absorber through which the particle passed.
  */
float CLoss::getEin(float energy, float thick)
{
  float dthick = 0.01;
  float de;
  float Ein= energy;

  //take small steps of thickness to increase the energy of fragment
  while(thick > 0.0){
    float thickness = min(thick,dthick);
    de = getDedx(Ein);
    Ein += de*thickness;
    if (Ein < 0. || Ein != Ein) return -1;
    thick -= dthick;
  }

  /* old loop
  for(;;){
    float thickness = min(thick,dthick);
    de = getDedx(Einput);
    Einput += de*thickness;
    if (thickness == thick) break;
    thick -= dthick;
  }
  */
  return Ein;
}


