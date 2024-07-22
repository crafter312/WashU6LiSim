#include <stdexcept>
#include <sstream>
#include <iomanip>

void plot_efficiencies(){
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
  
  TGraph *gr[Ex_num];// = new TGraph();
  TLegend *legend = new TLegend(0.15,0.7, 0.31, 0.9);
  
  int index = 0;
  double maxEff[Ex_num];
  double bestdist[Ex_num];

  for (int Ex_i=0; Ex_i<Ex_num; Ex_i++){
    gr[Ex_i] = new TGraph();
    maxEff[Ex_i] = 0;
    bestdist[Ex_i] = 0;
    
    for (int dist_i=0; dist_i<dist_num; dist_i++){
      index = Ex_i*dist_num + dist_i;
      tree->GetEntry(index);
      cout << index << " " << eff << " " << Ex << " " << distanceFromTarget << endl;
      gr[Ex_i]->SetPoint(dist_i, distanceFromTarget, eff);
      
      //find max eff
      if (eff > maxEff[Ex_i]){
        maxEff[Ex_i] = eff;
        bestdist[Ex_i] = distanceFromTarget;
      }
    }
  }
  
  for (int i=0; i<Ex_num; i++){
    cout << "Seperation Energy of 10." << i <<  " with max Efficiency of " << maxEff[i] << " at " << bestdist[i] << "mm" << endl;;
  }
  
  gr[0]->SetTitle("Efficiency of n+(He6+p)");
  gr[0]->GetXaxis()->SetTitle("Detector Distance (mm)");
  gr[0]->GetYaxis()->SetTitle("Efficiency");
  gr[0]->GetXaxis()->SetTitleSize(0.05);
  gr[0]->GetYaxis()->SetTitleSize(0.05);
  gr[0]->GetXaxis()->SetLabelSize(0.05);
  gr[0]->GetYaxis()->SetLabelSize(0.05);
  gr[0]->GetXaxis()->CenterTitle();
  gr[0]->GetYaxis()->CenterTitle();
  gr[0]->GetYaxis()->SetRangeUser(0,0.3);
  
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
  
  for (int Ex_i=0; Ex_i<Ex_num; Ex_i++){
    ostringstream leg;
    leg << "Ex = 10." << Ex_i;
    legend->AddEntry(gr[Ex_i], leg.str().c_str());
  }
  legend->SetTextSize(0.04);
  
  legend->Draw();
  
  
  for (int i=0; i<Ex_num; i++){
  
    ostringstream msg;
    msg << fixed << setprecision(3) << "#scale[0.8]{Ex = 10." << i << " effMax = " << maxEff[i] << " @ " << (int)bestdist[i] << "mm}";
    //TLatex *txt = new TLatex(170, 0.0055-0.0012*i, msg.str().c_str());
    TLatex *txt = new TLatex(170, 0.11-0.02*i, msg.str().c_str());
    txt->Draw();
  }

  
  //c1->SaveAs("plotting/efficienciesmulti_neut.gif");
  
  if (save) c1->Print("plotting/efficienciesmulti.png", "png");
  if (save) c1->Print("plotting/efficienciesmulti.svg", "svg");
}


