#include "fit.h"
#include <cmath>
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TArrow.h"
fit::fit(int n0, int NN0,int Npeaks0, double *para20, bool *useit0):minimizeND(n0)
{
  Npeaks = Npeaks0;

  Peak = new peak[Npeaks];

  para2 = para20;
  useit = useit0;
  n = n0;
  NN = NN0;
  const int vel = 2;  //2=all

  //read in efficiency data from a root file
  readeff();

  TFile fileX("/home/Li7star/Li7/sort_Li_CD2.root");
  TH1I* hist;

  hist = (TH1I*) fileX.Get("InvMass/7Li/Ex_7Li_ta");
  TH1I* hist_og = (TH1I*)hist->Clone();

  effcorr(hist);

  xlow = 2.5;
  xhigh = 6.0;
  low = hist->FindBin(xlow);
  high = hist->FindBin(xhigh)+1;
  nn = 0;
  ymax = 0.;
  for (int j=low;j<high;j++)
  {
    int i = j - low;
    yexp[i] = max(hist->GetBinContent(j),0.0);
    if (yexp[i] > ymax) ymax = yexp[i];
    x[i] = hist->GetBinCenter(j);
    syexp[i] = sqrt(max(hist_og->GetBinContent(j),0.0))*eff[j];
    if (yexp[i]<=0) syexp[i] = 1.;
    //syexp[i] = 1.;
    nn++;
    //if (x[i] > 14) break;
  }
  //take off the last point or so to make plot look better
  nn--;


  TFile * fsim1;
  fsim1 = new TFile ("/home/Li7star/li7sim/rootout/sim_ta_4642_090.root");
  TH1I* sim1 = (TH1I*) fsim1->Get("Ex");
  effcorr(sim1);

  //TFile * fsim2;
  //fsim2 = new TFile ("/home/Li7star/li7sim/rootout/sim_Li7_7500.root");
  //TH1I* sim2 = (TH1I*) fsim2->Get("Ex");

  for (int j=low;j<high;j++)
  {
    int i = j - low;
    if (i == nn) break;
    Peak[0].y[i] = sim1->GetBinContent(j);
    cout << "i" << i << " Peak[0].y[i]" << Peak[0].y[i] << endl;
    //Peak[1].y[i] = fake->GetBinContent(j);
    //Peak[1].y[i] = sim2->GetBinContent(j);
  }

  fsim1->Close();

  count = 0;
}

//overloads the funciton in the baseclass
//used by linmin in the base class minimizeND.
//Provides a 1D function for use by the base to call and get value of function
double fit::functND(double *para)
{
  int ii = 0;
  for (int  i = 0;i<NN;i++)
  {
    if (useit[i])
    {  
      para2[i] = para[ii];
      ii++;
    }
    //cout << para2[i] << " " << endl;
  }

  //cout << endl;
  //if (count == 10) abort();
  count++;

  double total = 0.;
  double totalCountsX = 0.;
  double totalFit = 0.;
  for (int i=0;i<nn;i++)
  {
    double k = getValue(i,para2);
    double delta = yexp[i]-k;
    if (yexp[i] == 0) continue;
    total += pow(delta/syexp[i],2);
  }



  if (std::isnan(total))
  {
    cout << " nan in fit" << endl;
    for (int i=0;i<NN;i++) cout << i << " " << para[i] << endl;
    abort();
  }

  return total;
}

void fit::read1new(string fn1)
{
  filename1 = fn1;
  cout << "file1: " << filename1 << endl;
  
  TFile * fsim1;
  fsim1 = new TFile (filename1.c_str());
  //fsim1 = new TFile ("/home/Li7star/Li7/plotting/out_sim_11300_175.root");
  TH1I* sim1 = (TH1I*) fsim1->Get("Ex");
  effcorr(sim1);


  for (int j=low;j<high;j++)
  {
    int i = j - low;
    if (i == nn) break;
    Peak[0].y[i] = sim1->GetBinContent(j);
  }

  fsim1->Close();
  count = 0;
}


//***********************************************
double fit::getValue(int i,double *para3)
{
  double back = getBackground(i,para3);

  double value =  0.;
  for (int k=0;k<Npeaks;k++) value += pow(para3[k],2)*Peak[k].y[i]*1.e-6;
      
  value = value + back;
  if (back < 0.) value += 200.;
  return value;
}

//readeff and effcorr are both used to accout for the efficiency of the Gobbi setup
void fit::readeff()
{
  //Get the efficiency information from a uniform Ex simulation
  TFile file_bkg("/home/Li7star/li7sim/rootout/sim_ta_uniform.root");
  TH1I* hist_eff = (TH1I*) file_bkg.Get("Ex");

  double max = 0;
  //find max value in response to uniform energy
  for (int i=0;i<400;i++)
  {
    if (hist_eff->GetBinContent(i) > max) max = hist_eff->GetBinContent(i);
  }
  //normalize efficiency to max value
  for (int i=0;i<400;i++)
  {
    if (hist_eff->GetBinContent(i) > 0)
    {
      eff[i] = max/hist_eff->GetBinContent(i);
    }
    else
    {
      eff[i] = 0;
    }
  }
}

void fit::effcorr(TH1I* h)
{
  //set bin content to account for efficiency
  for (int i=0;i<400;i++)
  {
    h->SetBinContent(i, (float)h->GetBinContent(i)*eff[i]); 
    cout << "effcorr: " << (float)h->GetBinContent(i) << " * " << eff[i] << "  = " << h->GetBinContent(i) << endl;
  }

}



//******************************************
double fit::getBackground(int i, double *para3)
{
  //return para3[14]*exp(-pow((x[i]-para3[15])/(para3[16]/10.),2)/2.);
  //return para3[16]*exp(-pow((x[i]-para3[17])/para3[18],2)/2.);
  int nn = Npeaks;

  //cout << "in getBackground" << para3[nn] << ", " << para3[nn+1] << ", " << para3[nn+2] << ", " << para3[nn+3] << endl;
  //abort();


  return (para3[nn]+para3[nn+3]*(x[i]-para3[nn+1])+ para3[nn+4]*pow(x[i]-para3[nn+1],2))/(1.+exp(-(x[i]-para3[nn+1])/para3[nn+2]));
  //return (para3[nn]+para3[nn+3]*(x[i]-para3[nn+1]) + para3[nn+4]*pow(x[i]-para3[nn+1],2) + para3[nn+5]*pow(x[i]-para3[nn+1],3))/(1.+exp(-(x[i]-para3[nn+1])/para3[nn+2]));

}

//********************************************
void fit::plot(double *para)
{
  ofstream  fileOut;
  TFile* fout;
  fout = new TFile("fit.root","RECREATE");
  fileOut.open("out.dat");

  TStyle * Sty = new TStyle("MyStyle","MyStyle");
  Sty->SetOptTitle(0);
  Sty->SetOptStat(0);
  //Sty->SetPalette(8,0);
  Sty->SetCanvasColor(10);
  Sty->SetCanvasBorderMode(0);
  Sty->SetFrameLineWidth(3);
  Sty->SetFrameFillColor(10);
  Sty->SetPadColor(10);
  Sty->SetPadTickX(1);
  Sty->SetPadTickY(1);
  Sty->SetPadBottomMargin(.18);
  Sty->SetPadLeftMargin(.18);
  Sty->SetPadTopMargin(.05);
  Sty->SetPadRightMargin(.03);
  Sty->SetHistLineWidth(3);
  Sty->SetHistLineColor(kRed);
  Sty->SetFuncWidth(3);
  Sty->SetFuncColor(kGreen);
  Sty->SetLineWidth(3);
  Sty->SetLabelSize(0.06,"xyz");
  Sty->SetLabelOffset(0.02,"y");
  Sty->SetLabelOffset(0.02,"x");
  Sty->SetLabelColor(kBlack,"xyz");
  Sty->SetTitleSize(0.08,"yz");
  Sty->SetTitleSize(0.06,"x");
   Sty->SetTitleOffset(1.16,"y");
  Sty->SetTitleOffset(1.2,"x");
  Sty->SetTitleFillColor(10);
  Sty->SetTitleTextColor(kBlack);
  Sty->SetTickLength(.05,"xz");
  Sty->SetTickLength(.025,"y");
  Sty->SetNdivisions(5,"xyz");
   Sty->SetLineWidth(2);
  Sty->SetEndErrorSize(0);
  gROOT->SetStyle("MyStyle");
  gROOT->ForceStyle();

  TCanvas canvas("fit");
  // canvas.SetLogy();
  gStyle->SetLineWidth(4);

  TH2S frame("frame","",10,xlow,xhigh,10,0,ymax*1.1);
  //TH2S frame("frame","",10,1.9,2.6,10,0,ymax*1.1);
  frame.GetXaxis()->SetTitle("E_{T} (MeV)");
  frame.GetXaxis()->SetNdivisions(410);
  frame.GetYaxis()->SetTitle("Counts");
  frame.GetXaxis()->CenterTitle();
  //frame.GetXaxis()->SetNdivisions(204,kTRUE);
  frame.GetYaxis()->CenterTitle();
  frame.Draw();

  double sx[500]={0.};
  TGraphErrors  gexp(nn,x,yexp,sx,syexp);
  gexp.SetMarkerStyle(20);
  gexp.SetMarkerSize(1.2);
  gexp.SetLineWidth(1.1);
  gexp.Draw("PE");

  double ytot[500];

  double maxy = 0.;
  int maxi = 0;
  for (int i=0;i<nn;i++)
  {
    ytot[i] = getValue(i,para);
    if (ytot[i] > maxy)
    {
      maxy = ytot[i];
      maxi = i;
    }
  }
  
  cout << " maxi = " << maxi << endl;


  double yb[500];
  for (int i=0;i<nn;i++)
  {
    yb[i] = getBackground(i,para);
  }
  TGraph  gb(nn,x,yb);
  gb.SetLineColor(4);
  gb.SetLineStyle(9);
  gb.Draw("C");


  double y1_s[500];
  double y2_s[500];
  double y3_s[500];
  double y4_s[500];


  for (int k=0;k<Npeaks;k++)
  {
    double sumy = 0.;
    for (int i=0;i<nn;i++)
    {
      Peak[k].ys[i] = Peak[k].y[i]*pow(para[k],2)*1.e-6;
      sumy += Peak[k].ys[i];
    }
    cout << " sum for peak " << k << " = " << sumy << endl;
  }

  TGraph * g[Npeaks];

  for (int k=0;k<Npeaks;k++)
  {
    g[k] = new TGraph(nn,x,Peak[k].ys);
    g[k]->SetLineColor(1);
    g[k]->SetLineWidth(2);
    g[k]->SetLineStyle(2);
    g[k]->Draw("L");
  }

  TGraph  gtot(nn,x,ytot);
  gtot.SetLineWidth(3);
  gtot.SetLineColor(2);
  gtot.Draw("L");


  for (int i=0;i<nn;i++)
  {
    fileOut <<  x[i] << " " <<  yexp[i] << " " << syexp[i] << " " <<  ytot[i] << " " <<  yb[i] << " " ;
    for (int k=0;k<Npeaks;k++) fileOut << Peak[k].ys[i] << " ";
    fileOut << endl;
  }


  canvas.Write();
  fout->Write();
}


//********************************
fit::~fit()
{
  delete Peak;
}
