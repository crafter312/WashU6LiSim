#ifndef _frag
#define _frag

#include "loss.h"
#include <string>
#include <memory>
#include "polyScat.h"
#include "random.h"
#include "frame.h"
//#include "fragment.h"
#include "Gobbiarray.h"
#include "tele.h"
using namespace std;

/**
 *!\brief information about a single fragment and its interection with the detector
 *
 */
class CFrag
{
 public:
  static CRandom ran;
  static float const pi; //!< 3.14159....

  float Z;
  float mass;
  float CsI_res;   

	//(float Z0, float mass0, string lossfile_CD2, string lossfile_Si, float CsI_res0, float thickness, Gobbiarray* array0, float scaleSmallAngle0, bool einstein0, bool useRealP0)
  CFrag(float,float,string,string,float,float,shared_ptr<Gobbiarray>,float=1.,bool=false,bool=false);
  ~CFrag();
  shared_ptr<Gobbiarray> Array;

  int hit(float,float,float);
  void getStripHit(int*, int*, int);
  int is_hit;

  void AddVelocity(double*);
  float Eloss(float); // energy loss in target
  float Egain(float); // corrects for energy loss in target
  void MultiScat(float);
  //void findVectors(float*,float*,float*,float*,float*,float,float,float,float,float,float,float,float,float);

  bool targetInteraction(float,float);
  bool SiliconInteraction();
  bool stopped;
  double FrontEnergy;
  double DeltaEnergy;

  CLoss *loss_C;
  CLoss *loss_Si;
  polyScat *pScat;

  CFrame *real;      //<!real particles energy, velocity, etc
  CFrame *recon;    //<!reconstructed properties


  float scaleSmallAngle;


	bool einstein;
  bool useRealP;   //<! use the real momentum, instead of detected value
  bool extra;

 private:
  float EgainHelper(float, CFrame*);
  
};

#endif
