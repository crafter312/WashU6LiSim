#include <iostream>
#include "Gobbiarray.h"
#include "frag.h"
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
  double Ebeam, Ex, gamma;

  //check if arguments were given on command line, if not defualt values are given
  if (argc > 1){Ebeam = atof(argv[1]);}
  else {Ebeam = 38.6;} //MeV
  if (argc > 2){Ex = atof(argv[2]);}
  else {Ex = 11.24;} //MeV
  if (argc > 3){gamma = atof(argv[3]);}
  else {gamma = 0.260;} //MeV
  cout << "argv[1]=E=" << Ebeam << "  argv[2]=Ex=" << Ex << "  argv[3]=gamma=" << gamma << endl;

  double distanceFromTarget = 239.5; //mm
  double Q =  Mass_6He + Mass_p - Mass_7Li;
  cout << "Q " << Q << endl;

  bool einstein = 1; //switch for newtonian(0) or relativistic(1) kinematics

  float thickness = 3.2; //mg/cm2  //target thickness
  float CsiRes = .01;    //resolution of Csi not needed for Si-Si;

  float const targetSize = 10.0;  //diameter of beam spot size in mm
  int Nevents = 50000;  // events to simulation

  bool useRealP = false; //true means use real angle and energies of fragment
                         //event reconstruction, to check effect of
                         //detector resolution

  //file for energy loss
  string Loss_H_in_CD2("Hydrogen_CD2.loss");
  string Loss_He_in_CD2("Helium_CD2.loss");

  //TODO need to update these so that they are in Silicon, currently a place holder
  string Loss_H_in_Si("Hydrogen_Si.loss");
  string Loss_He_in_Si("Helium_Si.loss");


  const int Nfrag = 2; //number of decay fragments
  CFrag ** frag;
  frag = new CFrag *[Nfrag];


  //scales the magnitude of small angle scattering
  float scale = 1;

  frag[0] = new CFrag(1., Mass_p/m0, Loss_H_in_CD2, Loss_H_in_Si, CsiRes, thickness, distanceFromTarget, scale, useRealP);
  frag[1] = new CFrag(2., Mass_6He/m0, Loss_He_in_CD2, Loss_He_in_Si, CsiRes, thickness, distanceFromTarget, scale, useRealP);

  CFrag *fragBeam;
  fragBeam = new CFrag(2., Mass_6He/m0, Loss_He_in_CD2, Loss_He_in_Si, CsiRes, thickness, distanceFromTarget, scale, useRealP);
  //create fragment to keep track of neutron fragment
  //CFragneut *fragneut;
  //fragneut = new CFragneut(Mass_n/m0);

  CDecay decay(Nfrag,frag,einstein);

  //define file for root spectra
  TFile * file = new TFile("sim.root","RECREATE");

  //define spectra
  //primary distribution of parent
  TH1F * hist_vel_P = new TH1F("vel_P","vel",100,1.5,4);
  TH1F * hist_theta_P = new TH1F("theta_P","theta",200,0,30);
  TH1F * hist_phi_P = new TH1F("phi_P","phi",180,0,360);
  //TH1F * hist_theta_neut_P = new TH1F("theta_neut_P","theta_neut",180,0,180);
  TH2F * kinematic_circle = new TH2F("kinematic_circle", "kin_circle", 200, -1, 1, 120, 2, 3.2);
  TH1F * hist_theta_beam_P = new TH1F("hist_theta_beam_P","theta",200,0,30);

  //reconsrtucted seconday distributions of parent;
  TH1F * hist_vel_S = new TH1F("vel_S","vel",100,1.5,4);
  TH1F * hist_theta_S= new TH1F("theta_S","theta",200,0,30);
  TH1F * hist_phi_S= new TH1F("phi_S","phi",180,0,360);
  //TH1F * hist_theta_neut_S = new TH1F("theta_neut_S","theta_neut",180,0,180);
  TH1F * hist_theta_beam_S_sharp = new TH1F("hist_theta_beam_S_sharp","theta",200,0,30);
  TH1F * hist_theta_beam_S_recon = new TH1F("hist_theta_beam_S_recon","theta",200,0,30);

  TH2I * DEE = new TH2I("DEE","",500,0,80,800,0,22); //E is x, DE is y
  TH1F * protonenergy = new TH1F("protonenergy","", 100,0,10);
  TH1F * He6energy = new TH1F("He6energy","", 100,0,40);

  //reconsrtucted Tertiary distributions of parent;
  //TH1F * hist_vel_T = new TH1F("vel_T","vel",100,1.5,4);
  //TH1F * hist_theta_T = new TH1F("theta_T","theta",200,0,30);
  //TH1F * hist_phi_T = new TH1F("phi_T","phi",180,0,360);
  //TH1F * hist_theta_neut_T = new TH1F("theta_neut_T","theta_neut",180,0,180);

  const double Ex_min = 9.5;
  const double dEx = 3.0;

  //TH1F * hist_neut_E_P = new TH1F("neut_E_P","",200, 0,20);
  //TH1F * hist_neut_E_S = new TH1F("neut_E_S","",200, 0,20);
  //TH1F * hist_neut_E_T = new TH1F("neut_E_T","",200, 0,20);

  TH1F * cos_thetaH = new TH1F("cos_thetaH","",100,-1,1);
  TH2F * hist_Erel_thetaH = new TH2F("Erel_thetaH","",200,0,8,25,-1,1);

  TH1I * hist_Erel_P = new TH1I("hist_Erel_P","",200,0,8);
  TH1I * hist_Erel = new TH1I("Erel","",200,0,8);
  TH1I * hist_Ex = new TH1I("Ex","",400,10,18);
  TH1I * hist_Ex_trans = new TH1I("Ex_trans","",400,10,18);
  TH1I * hist_Ex_trans_narrow = new TH1I("Ex_trans_narrow","",400,10,18);

  //maps of x-y positions of detected fragments
  TH2S * protonXY_S = new TH2S("protonXY_S","protonxy",100,-10,10,100,-10,10);
  TH2S * coreXY_S = new TH2S("coreXY_S","alphaxy",100,-10,10,100,-10,10);
  
  int Nstuck = 0;
  int Ndet = 0;
  int Nbeamscat = 0;

  //initialize Correlations so that it can read in the cross section data
  string Xsecfile = "../Li7FrescoCalc/Li7out_saved/new_36_fort.202";
  string elastic_correction = "../Li7FrescoCalc/HeCElasticout/fort16_simple.16";
  string elasXsecfile = "He6elasticC12.txt";

  Correlations *sampler;
  //initiallizing the Correlations class reads in the CM cross section from the file
  //and uses that to select a randomized value for phi and theta
  sampler = new Correlations(Xsecfile, elasXsecfile, elastic_correction, Ebeam, Ex);


  double pc0 = sqrt(pow(Ebeam+Mass_6He,2) - pow(Mass_6He,2));
  //cout << pc0 << endl;
  double P_acceptance = .012; //1.2% acceptence in MARS

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

    //add kinematics to beam
    //simulate MARS having +-1.2% acceptance range
    double pc = pc0*(1.-P_acceptance/2.) + decay.ran.Rndm()*P_acceptance*pc0;
    double Ebeam = sqrt(pow(pc,2) + pow(Mass_6He,2)) - Mass_6He;

    fragBeam->real->theta = sampler->thetaElastic;
    fragBeam->real->phi = sampler->phi;
    fragBeam->real->energy = Ebeam; //~6MeV/u He-6
    fragBeam->real->v[0] = sin(sampler->thetaElastic)*cos(sampler->phi);
    fragBeam->real->v[1] = sin(sampler->thetaElastic)*sin(sampler->phi);
    fragBeam->real->v[2] = cos(sampler->thetaElastic);

    //save info on beam
    hist_theta_beam_P->Fill(sampler->thetaElastic*180/acos(-1.));

    //determine if the beam hits the detector
    fragBeam->targetInteraction(dthick,thickness);
    fragBeam->SiliconInteraction();

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
      //loss of energy for each fragment through Silicon detector
      frag[i]->SiliconInteraction();
    }

    if (frag[0]->FrontEnergy > 15)
    {
      frag[0]->real->energy = 10 + (5*decay.ran.Rndm());
      frag[0]->FrontEnergy = frag[0]->real->energy;
    }


    //detect fragments, and skip if not all fragments detected
    int nhit = 0;
    int ishit = 0;

    for (int i=0;i<Nfrag;i++)
    {
      ishit = frag[i]->hit(xTarget,yTarget);
      nhit += ishit;
      if (ishit)
      {
        DEE->Fill(frag[i]->FrontEnergy,frag[i]->DeltaEnergy);
      }
      if (ishit == -1)
      {
        Nstuck++;
      }
    }

    if (nhit != Nfrag) continue;
    //cout << frag[0]->recon->energy << endl;

    //taken out but not tested, please check
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

    hist_Erel->Fill(Erel_S);
    hist_Ex->Fill(Ex_S);
    hist_Erel_thetaH->Fill(Erel_S,decay.cos_thetaH);
    cos_thetaH->Fill(decay.cos_thetaH);

    //look at transverse emiison for better resolutions

    if (fabs(decay.cos_thetaH) < 0.7) hist_Ex_trans->Fill(Ex_S);
    if (fabs(decay.cos_thetaH) < 0.5) hist_Ex_trans_narrow->Fill(Ex_S);


    hist_vel_S->Fill(decay.plfRecon->velocity);
    hist_theta_S->Fill(decay.plfRecon->theta*180./acos(-1.));
    hist_phi_S->Fill(decay.plfRecon->phi*180./acos(-1.));

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
  //delete fragneut;

  //beep at me when finished (sadly doesn't work anymore)
  //cout << "\a";
}
