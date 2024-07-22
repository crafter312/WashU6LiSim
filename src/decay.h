#include "frag.h"
#include "random.h"
#include "constants.h"
#include "profile.h"

#include <valarray>

/**
 *!\brief selects the veloity vectors of the secondary fragments
 */

class CDecay
{
  public:
    static CRandom ran;
    bool einstein;

    float ET;
    float cos_thetaH;

    CFrame *real[5]; 
    CFrame *recon[5];
    CFrame *plfRecon;
    CFrame *plfRecon2;
    CFrame *partCM[5];
    CFrag  **frag;  //!< information about the decay fragments
    profile *prof;

    int Nfrag;

    CDecay(int,CFrag**,bool einstein0);
    void GenerateProfile(double Ex, double Q);
    ~CDecay();
    float getErelReal();
    float getErelRecon();
    float getErel(CFrame**);
    float getErelNewton(CFrame**);
    float getErelRel(CFrame**);
    float getErel_at(CFrame*,CFrame*);

    void Mode2Body(double Ex, double gamma, double Q);
    void Mode2BodyExact(double Ex, double gamma, double Q);
		void ModeMicroCanonical(double Ex, double gamma, double Q);
    void ModeLineShapes();

    float sumA;
    float ErelRecon; //!<reconstructed relative kinetic energy
    
    double mass1;
    double mass2;

    int protonbranch;
    int neutronbranch;

};
