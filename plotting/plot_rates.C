#include <stdexcept>
#include <sstream>
#include <iomanip>

void plot_rates(){
  //plot options
  bool plot_neut = false;
  bool save = true;
  
  TCanvas *c1 = new TCanvas();
  c1->SetWindowSize(1600,1200);
  c1->SetRightMargin(0.09);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  gPad->SetTickx();
  gPad->SetTicky();
  
  TFile *input = new TFile("/home/n.dronchi/Documents/li7sim/simmulti_dist.root", "read");
  
  TTree *tree = (TTree*)input->Get("tree");
  
  double eff,eff_neut,Ex,distanceFromTarget;
  int entries = tree->GetEntries();
  
  tree->SetBranchAddress("eff", &eff);
  tree->SetBranchAddress("eff_neut", &eff_neut);
  tree->SetBranchAddress("Ex", &Ex);
  tree->SetBranchAddress("distanceFromTarget", &distanceFromTarget);
  
  const int Ex_num = 4;
  const int dist_num = 60;
  
  if (Ex_num*dist_num != entries) throw invalid_argument("Ex_num*dist_num != entries");
  
  
  int plot_num;
  float xs,xf,ys,yf;
  if (plot_neut){
    plot_num = Ex_num*2;
    xs = 0.15; ys = 0.65; xf = 0.4; yf = 0.9;
  }
  else{
    plot_num = Ex_num;
    xs = 0.15; ys = 0.7; xf = 0.31; yf = 0.9;
  }
  TGraph *gr[plot_num];
  TLegend *legend = new TLegend(xs,ys,xf,yf);
  
  //everything required to calculate reaction rate
  double density = 1.0; // g/cm^3
  double flux = 20000; // 1/s
  double Areal = 0.0006; // g/cm^2
  double mmass = 16.0; // g/mol
  double NA = 6.02; // *10^23 nuclei/mol
  double sigma = 104.57; //mb

  double power = pow(10,-4); // 10^23 * (1mb = 10^-27 cm^2)
  double rate; //hr-1
  
  
  int index = 0;

  for (int Ex_i=0; Ex_i<Ex_num; Ex_i++){
    gr[Ex_i] = new TGraph();
    
    for (int dist_i=0; dist_i<dist_num; dist_i++){
      index = Ex_i*dist_num + dist_i;
      tree->GetEntry(index);
      //calc rate accounting for efficiency
      rate = flux*2*density*NA/mmass*sigma*Areal/density*power*3600*eff;
      //cout << index << " " << rate << " " << Ex << " " << distanceFromTarget << endl;
      gr[Ex_i]->SetPoint(dist_i, distanceFromTarget, rate);
    }
  }
  //optional plot for neutron efficiency
  if (plot_neut){
    for (int Ex_i=0; Ex_i<Ex_num; Ex_i++){
      gr[Ex_i+Ex_num] = new TGraph();

      for (int dist_i=0; dist_i<dist_num; dist_i++){
        index = Ex_i*dist_num + dist_i;
        tree->GetEntry(index);
        //calc rate accounting for efficiency
        rate = flux*2*density*NA/mmass*sigma*Areal/density*power*3600*eff_neut;
        //cout << index << " " << rate << " " << Ex << " " << distanceFromTarget << endl;
        gr[Ex_i+Ex_num]->SetPoint(dist_i, distanceFromTarget, rate);
      }
    }
  }
  
  gr[0]->SetTitle("Rate of n+(He6+p)");
  gr[0]->GetXaxis()->SetTitle("Detector Distance (mm)");
  gr[0]->GetYaxis()->SetTitle("Rate (hr^{-1})");
  gr[0]->GetXaxis()->SetTitleSize(0.05);
  gr[0]->GetYaxis()->SetTitleSize(0.05);
  gr[0]->GetXaxis()->SetLabelSize(0.05);
  gr[0]->GetYaxis()->SetLabelSize(0.05);
  gr[0]->GetXaxis()->CenterTitle();
  gr[0]->GetYaxis()->CenterTitle();
  //gr[0]->GetYaxis()->SetRangeUser(0,0.02);
  
  gr[0]->SetLineColor(1);
  gr[0]->SetLineStyle(2);
  gr[0]->SetMarkerStyle(20);
  gr[0]->SetMarkerColor(1);
  gr[0]->SetMarkerSize(2);
  
  gr[0]->Draw("ACP");
  
  for (int Ex_i=1; Ex_i<Ex_num; Ex_i++){
    gr[Ex_i]->SetLineColor(Ex_i+1);
    gr[Ex_i]->SetMarkerColor(Ex_i+1);
    gr[Ex_i]->SetMarkerStyle(Ex_i+20);
    gr[Ex_i]->SetMarkerSize(2);
    gr[Ex_i]->SetLineStyle(2);
    gr[Ex_i]->Draw("CP");
  }
  if (plot_neut){
    for (int Ex_i=0; Ex_i<Ex_num; Ex_i++){
      gr[Ex_i+Ex_num]->SetLineColor(Ex_i+1);
      gr[Ex_i+Ex_num]->SetMarkerColor(Ex_i+1);
      gr[Ex_i+Ex_num]->SetMarkerStyle(Ex_i+24);
      if (Ex_i == 3) gr[Ex_i+Ex_num]->SetMarkerStyle(32);
      gr[Ex_i+Ex_num]->SetMarkerSize(1);
      gr[Ex_i+Ex_num]->SetLineStyle(2);
      gr[Ex_i+Ex_num]->Draw("CP");
    }
  }
  
  
  for (int Ex_i=0; Ex_i<Ex_num; Ex_i++){
    ostringstream leg;
    leg << "#scale[0.8]{Ex = 10." << Ex_i << "}";
    legend->AddEntry(gr[Ex_i], leg.str().c_str());
  }
  if (plot_neut){
    for (int Ex_i=0; Ex_i<Ex_num; Ex_i++){
      ostringstream leg;
      leg << "#scale[0.8]{Ex = 10." << Ex_i << "w/ Neut}";
      legend->AddEntry(gr[Ex_i+Ex_num], leg.str().c_str());
    }
  }
  
  legend->SetTextSize(0.04);
  legend->Draw();
  
  if (save && plot_neut) c1->Print("plotting/ratesmutli_neut.png", "png");
  if (save && plot_neut) c1->Print("plotting/ratesmutli_neut.svg", "svg");
  if (save && !plot_neut) c1->Print("plotting/ratesmutli.png", "png");
  if (save && !plot_neut) c1->Print("plotting/ratesmutli.svg", "svg");
  
}


