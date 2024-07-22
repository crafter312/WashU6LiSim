#include <iostream>
#include "Gobbiarray.h"
#include "Neutarray.h"
#include "frag.h"
#include "fragneut.h"
#include "decay.h"
#include "constants.h"
#include "correlations.h"

#include "TH1F.h"
#include "TH2S.h"
#include "TFile.h"
#include "TF1.h"

using namespace std;


int main(int argc, char *argv[])
{
  //define excitation energy of state, with, and Q valeu for break up
  double E, Ex, distanceFromTarget;

  //check if arguments were given on command line, if not defualt values are given
  if (argc > 1){E = atof(argv[1]);}
  else {E = 38.6;} //MeV
  if (argc > 2){Ex = atof(argv[2]);}
  else {Ex = 11.24;} //MeV
  if (argc > 3){distanceFromTarget = atof(argv[3]);}
  else {distanceFromTarget = 241;} //mm
  cout << "argv[1]=E=" << E << "  argv[2]=Ex=" << Ex << "  argv[3]=distanceFromTarget=" << distanceFromTarget << endl;

  double gamma = 0.0;
  double Q =  Mass_6He + Mass_p - Mass_7Li;
  cout << "Q " << Q << endl;

  bool einstein = 1; //switch for newtonian(0) or relativistic(1) kinematics

  float thickness = 3.0; //mg/cm2  //target thickness
  float CsiRes = .05;    //resolution of Csi not needed for Si-Si;

  float const targetSize = 4.;  //diameter of beam spot size in mm
  int Nevents = 50000;  // events to simulation

  bool useRealP = false; //true means use real angle and energies of fragment
                         //event reconstruction, to check effect of
                         //detector resolution

  //file for energy loss
  string Loss1("Hydrogen_CD2.loss");
  string Loss2("Helium_CD2.loss");

  const int Nfrag = 2; //number of decay fragments
  CFrag ** frag;
  frag = new CFrag *[Nfrag];


  //scales the magnitude of small angle scattering
  float scale = 1;

  frag[0] = new CFrag(1.,Mass_p/m0,Loss1,CsiRes,thickness,distanceFromTarget,scale,useRealP);
  frag[1] = new CFrag(6.,Mass_6He/m0,Loss2,CsiRes,thickness,distanceFromTarget,scale,useRealP);

  CFrag *fragBeam;
  fragBeam = new CFrag(6.,Mass_6He/m0,Loss2,CsiRes,thickness,distanceFromTarget,scale,useRealP);
  //create fragment to keep track of neutron fragment
  CFragneut *fragneut;
  fragneut = new CFragneut(Mass_n/m0);

  CDecay decay(Nfrag,frag,einstein);

  //define file for root spectra
  TFile * file = new TFile("sim.root","RECREATE");

  //define spectra
  //primary distribution of parent
  TH1F * hist_vel_P = new TH1F("vel_P","vel",100,1.5,4);
  TH1F * hist_theta_P = new TH1F("theta_P","theta",200,0,30);
  TH1F * hist_phi_P = new TH1F("phi_P","phi",180,0,360);
  TH1F * hist_theta_neut_P = new TH1F("theta_neut_P","theta_neut",180,0,180);
  TH2F * kinematic_circle = new TH2F("kinematic_circle", "kin_circle", 200, -1, 1, 120, 2, 3.2);
  TH1F * hist_theta_beam_P = new TH1F("hist_theta_beam_P","theta",200,0,30);

  //reconsrtucted seconday distributions of parent;
  TH1F * hist_vel_S = new TH1F("vel_S","vel",100,1.5,4);
  TH1F * hist_theta_S= new TH1F("theta_S","theta",200,0,30);
  TH1F * hist_phi_S= new TH1F("phi_S","phi",180,0,360);
  TH1F * hist_theta_neut_S = new TH1F("theta_neut_S","theta_neut",180,0,180);
  TH1F * hist_theta_beam_S_sharp = new TH1F("hist_theta_beam_S_sharp","theta",200,0,30);
  TH1F * hist_theta_beam_S_recon = new TH1F("hist_theta_beam_S_recon","theta",200,0,30);

  TH1F * protonenergy = new TH1F("protonenergy","", 100,0,10);
  TH1F * He6energy = new TH1F("He6energy","", 100,0,40);

  //reconsrtucted Tertiary distributions of parent;
  TH1F * hist_vel_T = new TH1F("vel_T","vel",100,1.5,4);
  TH1F * hist_theta_T = new TH1F("theta_T","theta",200,0,30);
  TH1F * hist_phi_T = new TH1F("phi_T","phi",180,0,360);
  TH1F * hist_theta_neut_T = new TH1F("theta_neut_T","theta_neut",180,0,180);

  const double Ex_min = 9.5;
  const double dEx = 3.0;

  TH1F * hist_neut_E_P = new TH1F("neut_E_P","",200, 0,20);
  TH1F * hist_neut_E_S = new TH1F("neut_E_S","",200, 0,20);
  TH1F * hist_neut_E_T = new TH1F("neut_E_T","",200, 0,20);
  TH1F * hist_Ex_S = new TH1F("Ex_S","",400,10,18);
  TH1F * hist_Erel_S = new TH1F("Erel_S","",200,0,dEx);

  TH1F * cos_thetaH = new TH1F("cos_thetaH","",100,-1,1);
  TH2F * hist_Erel_thetaH = new TH2F("Erel_thetaH","",200,0,8,25,-1,1);

  TH1F * hist_Erel_P = new TH1F("Erel_P","",200,0,dEx);
  TH1F * hist_Ex_trans_S = new TH1F("Ex_trans_S","",200,Ex_min,Ex_min+dEx);
  TH1F * hist_Ex_trans_narrow_S = new TH1F("Ex_trans_narrow_S","",200,Ex_min,Ex_min+dEx);

  //maps of x-y positions of detected fragments
  TH2S * protonXY_S = new TH2S("protonXY_S","protonxy",100,-10,10,100,-10,10);
  TH2S * coreXY_S = new TH2S("coreXY_S","alphaxy",100,-10,10,100,-10,10);
  
  int Nstuck = 0;
  int Ndet = 0;
  int Ndet_neut = 0;
  int Nbeamscat = 0;

  //initialize Correlations so that it can read in the cross section data
  string Xsecfile = "../Li7FrescoCalc/Li7out_saved/new_36_fort.202";
  string elastic_correction = "../Li7FrescoCalc/HeCElasticout/fort16_simple.16";
  string elasXsecfile = "He6elasticC12.txt";

  Correlations *sampler;
  //initiallizing the Correlations class reads in the CM cross section from the file
  //and uses that to select a randomized value for phi and theta
  sampler = new Correlations(Xsecfile, elasXsecfile, elastic_correction, E, Ex);

  for (int index=0;index<Nevents;index++)
  {

    //progress updates
    //if (index%10000 == 0){cout << (float)index/(float)Nevents << endl;}
    if (index%1000 == 0) cout << '\xd' <<  index << " of " << Nevents << flush;

    // distance in target that produced has to pass to get out
    double dthick = thickness*decay.ran.Rndm();

    // beam spot at target
    double rTarget = sqrt(decay.ran.Rndm())*targetSize/2.;
    double theta = 2.*acos(-1.)*decay.ran.Rndm();

    float xTarget = rTarget*cos(theta);
    float yTarget = rTarget*sin(theta);


    //need to re-randomize the angles for each passthrough
    sampler->randomAngles();

    hist_vel_P->Fill(sampler->Vpplab);
    hist_theta_P->Fill((sampler->thetaLab)*180/acos(-1.));
    hist_phi_P->Fill(sampler->phi*180/acos(-1.));

    //add kinematics to neutron fragment class
    //double VVneut[3];
    //VVneut[0] = sampler->Vttlab*sin(sampler->thetaNeut)*cos(sampler->phineut);
    //VVneut[1] = sampler->Vttlab*sin(sampler->thetaNeut)*sin(sampler->phineut);
    //VVneut[2] = sampler->Vttlab*cos(sampler->thetaNeut);
    fragneut->setVelocity(VVneut);
    //save hist info on neutron
    hist_theta_neut_P->Fill(sampler->thetaNeut*180/acos(-1.));
    hist_neut_E_P->Fill(fragneut->neut->energy);

    //add kinematics to beam
    fragBeam->real->theta = sampler->thetaElastic;
    fragBeam->real->phi = sampler->phi;
    fragBeam->real->energy = 6*6; //~6MeV/u He-6
    fragBeam->real->v[0] = sin(sampler->thetaElastic)*cos(sampler->phi);
    fragBeam->real->v[1] = sin(sampler->thetaElastic)*sin(sampler->phi);
    fragBeam->real->v[2] = cos(sampler->thetaElastic);

    //save info on beam
    //cout << "thetaElastic " << sampler->thetaElastic << "  " <<sampler->thetaElastic*180/acos(-1.) << endl;
    hist_theta_beam_P->Fill(sampler->thetaElastic*180/acos(-1.));
    //hist_beam_vel_P->Fill(fragBeam->real->velocity);
    //cout << fragBeam->recon->velocity << endl;

    //determine if the beam hits the detector
    fragBeam->targetInteraction(dthick,thickness);
    int beamhit = fragBeam->hit(xTarget,yTarget);
    if (beamhit)
    {
      //hist_theta_beam_S->Fill(fragBeam->recon->theta*180/acos(-1.));
      hist_theta_beam_S_sharp->Fill(sampler->thetaElastic*180/acos(-1.));
      hist_theta_beam_S_recon->Fill(fragBeam->recon->theta*180/acos(-1.));
      Nbeamscat++;
    }

    
    //add kinematics to the charged particles
    //velocity vector, z axis is beam axis
    double VVparent[3];
    VVparent[0] = sampler->Vpplab*sin(sampler->thetaLab)*cos(sampler->phi);
    VVparent[1] = sampler->Vpplab*sin(sampler->thetaLab)*sin(sampler->phi);
    VVparent[2] = sampler->Vpplab*cos(sampler->thetaLab);

    kinematic_circle->Fill(VVparent[0], VVparent[2]);

    //decay parent fragment, add sets velocity vectors of fragments to the seperation
    decay.Mode2Body(Ex,gamma,Q);

    hist_Erel_P->Fill(decay.ET);

    //transfrom decay vectors to lab frame by adding initial velocity of partent Li7 to both fragments
    for (int i=0;i<Nfrag;i++){ frag[i]->AddVelocity(VVparent);}

    //interaction of fragements in target material. Calcs energy loss in target, change in scatter angle,
    //and wheter fragment is stopped within target
    for (int i=0;i<Nfrag;i++)
    {
      frag[i]->targetInteraction(dthick,thickness);
    }

    //detect fragments, and skip if not all fragments detected
    int nhit = 0;
    for (int i=0;i<Nfrag;i++)
    {
      nhit += frag[i]->hit(xTarget,yTarget);
      if (frag[i]->hit(xTarget,yTarget) == -1)
      {
        Nstuck++;
      }
    }
    if (nhit != Nfrag) continue;
    //cout << frag[0]->recon->energy << endl;
    if (frag[0]->recon->energy < 2.5) continue;

    //if seperation enery is small, make sure they hit different silicon strips
    int stripx[Nfrag];
    int stripy[Nfrag];
    //collect what strips are hit
    for (int i=0;i<Nfrag;i++)
    {
      frag[i]->getStripHit(stripx, stripy, i);
    }

    bool doublehit = false;
    //loop through all pairs of strips
    for (int i=0;i<Nfrag;i++)
    {
      for (int j=i+1; j<Nfrag; j++)
      {
        //check if it double hit strip
        if (stripx[i] == stripx[j] || stripy[i] == stripy[j])
        {
          doublehit = true; //use to skip after loop
        }
      }
    }

    if (doublehit) continue;
    
    //We have a detection
    Ndet++;

    for (int i=0; i<Nfrag; i++)
    {
      frag[i]->Egain(thickness/2.0);
    }


    protonenergy->Fill(frag[0]->recon->energy);
    He6energy->Fill(frag[1]->recon->energy);


    //get reconstructed  relative energy between fragements
    float Erel_S = decay.getErelRecon();

    //get reconstructed excitation energy

    float Ex_S = Erel_S+Q;

    hist_Ex_S->Fill(Ex_S);
    hist_Erel_S->Fill(Erel_S);
    hist_Erel_thetaH->Fill(Erel_S,decay.cos_thetaH);
    cos_thetaH->Fill(decay.cos_thetaH);

    //look at transverse emiison for better resolutions
    if (fabs(decay.cos_thetaH) < 0.5) hist_Ex_trans_S->Fill(Ex_S);
    if (fabs(decay.cos_thetaH) < 0.2) hist_Ex_trans_narrow_S->Fill(Ex_S);


    hist_vel_S->Fill(decay.plfRecon->velocity);
    hist_theta_S->Fill(decay.plfRecon->theta*180./acos(-1.));
    hist_phi_S->Fill(decay.plfRecon->phi*180./acos(-1.));

    //TODO modify so that you only get angle from the p-terp crystal
    //hist_theta_neut_S->Fill(sampler->thetaNeut*180/acos(-1));
    //hist_neut_E_S->Fill(fragneut->neut->energy);


    float x = frag[0]->recon->x/10.;
    float y = frag[0]->recon->y/10.;

    protonXY_S->Fill(x,y);

    x = frag[1]->recon->x/10.;
    y = frag[1]->recon->y/10.;
    coreXY_S->Fill(x,y);

  }

  cout << endl;
  cout << endl;

  cout << "(6He + p) coincidence efficiency = " << (float)Ndet/(float)Nevents << endl;
  cout << "n + (6He + p) coincidence efficiency = " << (float)Ndet_neut/(float)Nevents << endl;
  cout << "(6He + p) was stuck in target = " << Nstuck << " percent = " << (float)Nstuck/(float)Nevents << endl;

  cout << "6He elastic scatter det " << (float)Nbeamscat/(float)Nevents << endl;
  
  //make fit to measure resolution of Invarient mass peak
  //TF1 * fit = new TF1("fit", "gaus", Ex_min, Ex_min+dEx);
  //hist_Ex_S->Fit("fit", "R");
  //double mean = fit->GetParameter(1);
  //double sigma = fit->GetParameter(2);
  //cout << "mean " << mean << " sigma " << sigma << endl;


  //write out root file
  file->Write();

  //clean up, clean up
  for (int i=0;i<Nfrag;i++)
  {
    delete frag[i];
  }
  //everybody everywhere
  delete frag;
  //clean up, clean up
  delete sampler;
  //everybody do your share
  delete fragneut;

  //beep at me when finished
  cout << "\a";
}
