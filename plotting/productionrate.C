#include <vector>
#include <string>
#include "TGraph.h"

void get_jacobians(vector<double> &returnvec, double theta, double Vpp, double Vtt, double VCM){
    const double PI = 3.141592;
    // projections of projectile out on x and z in CM frame
    double Vxpp = Vpp*sin(theta);
    double Vzpp = Vpp*cos(theta);
    // projections of projectile out on x and z in Lab frame
    double Vxpplab = Vxpp;
    double Vzpplab = Vzpp + VCM ;
    // lab velocity and angle of the projectile
    double Vpplab = sqrt(pow(Vxpplab,2) + pow(Vzpplab,2));
    double labangle = acos(Vzpplab/Vpplab);
    // jacobian transformation for projectile fragment
    double Jacobian = pow(Vpplab/Vpp,2)/cos(labangle - theta);
    // target is always 180 deg away in center of mass frame
    double NangCM = PI - theta;
    // projections of projectile out on x and z in CM frame
    double Vxtt = Vtt*sin(PI - theta);
    double Vztt = Vtt*cos(PI - theta);
    // projections of projectile out on x and z in Lab frame
    double Vxttlab = Vxtt;
    double Vzttlab = Vztt + VCM;
    // lab velocity and angle of the projectile
    double Vttlab = sqrt(pow(Vxttlab,2) + pow(Vzttlab,2));
    double Nangle = acos(Vzttlab/Vttlab);
    // jacobian transformation for projectile fragment
    double Jacobian2 = pow(Vttlab/Vtt,2)/cos(Nangle - NangCM);
  
    returnvec.push_back(labangle);
    returnvec.push_back(Jacobian);
    returnvec.push_back(Nangle);
    returnvec.push_back(Jacobian2);
}    


void CMtoLab(vector<double> &th, vector<double> &Xsec, vector<double> &ang_Li, vector<double> &Xsec_Li,
             vector<double> &ang_n, vector<double> &Xsec_n, double E, double Ex){
    const double PI = 3.141592;
    double Mp = 6.02; // mass of projectile He6
    double Mt = 2.014; // mass of target H2
    double Mpp = 7.016; // mass of outgoing projectile Li7
    double Mtt = 1.008; // mass of outgoing target 
    double Mred = Mpp*Mtt/(Mpp+Mtt); // reduced mass

    double EnergyPA = E/Mp; // energy per nucleon
    
    double Vbeam = sqrt(2*EnergyPA)*0.983; // beam velocity
    
    double VCM = Vbeam * Mp/(Mp+Mt); // velocity of CM
    double VpCM = Vbeam-VCM; // velocity of projectile in CM frame
    double VtCM = VCM; // velocity of target in CM frame is velocity in CM
    
    double ECMin = 0.5*Mp*pow(VpCM*0.983,2) + 0.5*Mt*pow(VtCM/0.983,2);
    double Q = 7.75; // MeV
    double ECMout = ECMin + Q - Ex;
    
    // relative velocity between projectile and target in exit channel
    double Vrel = sqrt(ECMout/0.983*2/Mred);
    
    // center of mass velocity of projectile fragment in exit channel
    double Vpp = Vrel*Mtt/(Mtt+Mpp);
    // CM velocity of target fragment in exit channel
    double Vtt = Vrel - Vpp;
    
    double theta;
    vector<double> jac;
    // Loop over CM cross section to convert every point to Lab frame
    for (int i=0; i<th.size(); i++){
        // deg to radians
        theta = (th[i]*PI)/180;
        
        // get the jacobians and lab angles
        get_jacobians(jac, theta, Vpp, Vtt, VCM);
        // Lilabangle=jac[0], Lijacobian=jac[1],  nlabangle=jac[2], njacobian=jac[3]
        
        // save transformed coordinates
        ang_Li.push_back(jac[0]*180/PI);
        ang_n.push_back(jac[2]*180/PI);
        Xsec_Li.push_back(abs(Xsec[i]*jac[1]));
        Xsec_n.push_back(Xsec[i]*jac[3]);
        jac.clear();
    }
}

void readfile202(const string filename, vector<double> &th, vector<double> &Xsec){
    fstream file;
    file.open(filename, ios::in);
    
    string line;   
 
    // go through file and take out cross section data
    while (getline(file, line)){
        if (line.compare(0,1,"#") == 0){
            // print out the comments
            //cout << line << endl;
        }
        if (line.compare(0,3,"   ") == 0){   
            // store cross section data in vectors        
            th.push_back(stod(line.substr(3,5)));
            Xsec.push_back(stod(line.substr(15,5)));
        }
    }
    file.close();
}

void productionrate(){
  const double PI = acos(-1.);
  //TRandom2 *rand = new TRandom2(0);

  vector<double> th;
  vector<double> Xsec;
  vector<double> ang_Li;
  vector<double> Xsec_Li;
  vector<double> ang_n;
  vector<double> Xsec_n;

  //readfile202("Li7out_saved/gs_36_fort.202", th, Xsec);    
  //CMtoLab(th, Xsec, ang_Li, Xsec_Li, ang_n, Xsec_n, 36, 0);
  //TGraph *gr_gs_Li = new TGraph(ang_Li.size(), &ang_Li[0], &Xsec_Li[0]);
  //TGraph *gr_gs_n = new TGraph(ang_n.size(), &ang_n[0], &Xsec_n[0]);

  double density = 1.0; // g/cm^3
  double flux = 20000; // 1/s
  double Areal = 0.0006; // g/cm^2
  double mmass = 16.0; // g/mol
  double NA = 6.02; // *10^23 nuclei/mol

  double power = pow(10,-4); // 10^23 * (1mb = 10^-27 cm^2)

  readfile202("/home/n.dronchi/Documents/Li7Excited/Li7out_saved/new_36_fort.202", th, Xsec);
  CMtoLab(th, Xsec, ang_Li, Xsec_Li, ang_n, Xsec_n, 36, 10);
  
  //plot lab cross section
  TGraph *gr_new_Li = new TGraph(ang_Li.size(), &ang_Li[0], &Xsec_Li[0]);
  gr_new_Li->GetXaxis()->SetRangeUser(0,18);
  gr_new_Li->GetYaxis()->SetRangeUser(0,4000);
  auto *c1 = new TCanvas("c1", "Lab Li7 Differential Cross Sections");
  c1->SetGrid();
  gStyle->SetGridStyle(3);
  gr_new_Li->Draw("ACP");

  double sigma = 0; // temp value for sigma in millibarn
  double step;
  double theta;
  // 2PI is from azimuthal, integrate[ d(sigma_Li7)/dOmega *sin(theta)*dtheta] for all theta > 6deg
  for (int i=0; i<ang_Li.size()-1; i++){
    // step size changes in lab frame
    step = ((ang_Li[i+1] - ang_Li[i])*PI)/180.0;
    theta = (ang_Li[i]*PI)/180.0; //deg->rad
    sigma += 2*PI*Xsec_Li[i]*sin(theta)*abs(step);
  }
  
  cout << "total detected sigma " << sigma << "mb,   expected from fresco 104.57mb" << endl;
  //cout << " (Should be 512mb if full integration and gs)" << endl;


  double rate = flux*2*density*NA/mmass*sigma*Areal/density*power;
  cout << "we will detect " << rate << "per second, " << rate*60 << "per minute, or " << rate*60*60 << "per hour" << endl;
 
}


