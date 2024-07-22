/**
 *!\brief solves a system of linear equations
 */
class sle
{
public:
  float ** M; //!< points to a matrix containing the coefficents of the equations
  float * Y; //!< points to a vector containing the solution 
  int N; //!< number of equations
  sle(int);
  ~sle();
  void solve();
};

sle::sle(int N0)
{
  N = N0;
  M = new float*[N];
  for (int i=0;i<N;i++) M[i] = new float[N];
  Y = new float[N];
}
//**************************************************************************
  /*
   * destructor
   */
sle::~sle()
{
  delete [] Y;
  for (int i=0;i<N;i++) delete [] M[i];
  delete [] M;
}
//**************************************************************************
  /**
   * Solves the equation, after one manually enters the values of M and Y
   */
void sle::solve()
{
  for (int j=0;j<N;j++)
  {
    float constant = M[j][j];
    for (int i=j;i<N;i++) M[j][i] /= constant;
    Y[j] /= constant;
    for (int i=0;i<N;i++)
    {
      if (i != j)
      {
        constant = M[i][j];
        for (int k=j;k<N;k++) M[i][k] -= constant*M[j][k];
        Y[i] -= constant*Y[j];
      }
    }
  }
}


class chi2
{
public:
  int binlow;
  int binhigh;
  int Nvals;

  float *y;
  float *a;
  float *b;
  float *c;

  sle *LinEq;

  chi2(int, int, int, float*, float*);
  ~chi2();

  float CalcChisqr(float);
  //void fitparameters(float);
  //float CalcUncertainty(float, float, float);
};

chi2::~chi2()
{
  delete LinEq;
}


chi2::chi2(int Nvals0, int binlow0, int binhigh0, float *y0, float *a0)
{
  binlow = binlow0;
  binhigh = binhigh0;
  Nvals = Nvals0;
  y = y0; //data
  a = a0; //sim

  //Fitting only 1 function sim data.
  LinEq = new sle(1); 
}

float chi2::CalcChisqr(float fitA)
{
  //calc Starting chi^2
  float chi2val = 0;
  for (int i=binlow; i<binhigh; i++)
  {
    float sigma2 = y[i];
    if ( sigma2 == 0) { sigma2 = 1;}
    chi2val += pow((y[i] - fitA*a[i]),2)/sigma2;
  }
  //chi2val = chi2val/(Nvals-1); //divide by total number of values fit
  return chi2val;
}
/*
void chi2::fitparameters(float fixedB)
{
  //subtract out fixed value on the fitting
  float y_fixed[Nvals];
  for (int i=0; i<Nvals; i++)
  {
    y_fixed[i] = y[i] - fixedB*b[i];
  }

  //fill "sle" LinEq.Y[k], where y[i] is the solution.
  //k=0 is the bkg data
  float sum = 0;
  for (int i=binlow; i<binhigh+1; i++)
  {
    sum += y_fixed[i] * a[i];
  }
  LinEq->Y[0] = sum;
  //k=1 is the compton geant data
  sum = 0;
  for (int i=binlow; i<binhigh+1; i++)
  {
    sum += y_fixed[i] * c[i];
  }
  LinEq->Y[1] = sum;

  //start filling matrix M[Nfits][Nfits] = M[2][2]
//Get the diagonal elements
  sum = 0;
  for (int i=binlow; i<binhigh+1; i++)
  {
    sum += a[i] * a[i];
  }
  LinEq->M[0][0] = sum;
  sum = 0;
  for (int i=binlow; i<binhigh+1; i++)
  {
    sum += c[i] * c[i];
  }
  LinEq->M[1][1] = sum;
//get the off diagonal elements
  sum = 0;
  for (int i=binlow; i<binhigh+1; i++)
  {
    sum += a[i] * c[i];
  }
  LinEq->M[1][0] = sum;
  LinEq->M[0][1] = sum;

  LinEq->solve();
}


float chi2::CalcUncertainty(float fitA, float fitB, float fitC)
{
  float chi2 = CalcChisqr(fitA, fitB, fitC);

  float stepsize = 1./10000.*fitB;
  float curfitB = fitB;
  int step = 0; // if program isn't working check if this is too high
  float newchi2 = 0;
  //for different values of A
  while (newchi2 < chi2+1.0)
  {
    step++;
    float newfitA, newfitC;
    curfitB = fitB + stepsize*step;
    //cout << "curfitB " << curfitB << "   = " << fitB << " + " << stepsize*step << endl;

    //fit new values for B and C based on new A
    fitparameters(curfitB);

    newfitA = LinEq->Y[0];
    newfitC = LinEq->Y[1];
    
    newchi2 = CalcChisqr(newfitA,curfitB,newfitC);

    if (step > 10000) //protect against infinite loop
    {
      cout << "issue with CalcUncertainty!!!!" << endl;
      abort();
    }

    //cout << "fitA " << fitA << "," << newfitA << "   fitB " << fitB << "," << curfitB << "   fitC " << fitC << "," << newfitC << endl;
    //cout << "chi2 " << chi2 << " newchi2 " << newchi2 << endl;
  }

  if (step < 10) {cout << "low number of steps taken, think about decreasing stepsize" << endl;}

  return curfitB - fitB;
}
*/

void Fitspec()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  TStyle * Sty = (TStyle*)gROOT->FindObject("MyStyle");
  if(!Sty)
    {
      Sty = new TStyle("MyStyle","MyStyle");
    }
  Sty->SetOptTitle(0);
  Sty->SetOptStat(0);
  Sty->SetPalette(1,0);
  Sty->SetCanvasColor(10);
  Sty->SetCanvasBorderMode(0);
  Sty->SetFrameLineWidth(3);
  Sty->SetFrameFillColor(10);
  Sty->SetPadColor(10);
  Sty->SetPadTickX(1);
  Sty->SetPadTickY(1);
  Sty->SetPadBottomMargin(.17);
  Sty->SetPadTopMargin(.03);
  Sty->SetPadLeftMargin(.17);
  Sty->SetPadRightMargin(.05);
  Sty->SetHistLineWidth(3);
  Sty->SetHistLineColor(kBlack);
  Sty->SetFuncWidth(3);
  Sty->SetFuncColor(kRed);
  Sty->SetLineWidth(3);
  Sty->SetLabelSize(0.06,"xyz");
  Sty->SetLabelOffset(0.02,"y");
  Sty->SetLabelOffset(0.02,"x");
  Sty->SetLabelColor(kBlack,"xyz");
  Sty->SetTitleSize(0.06,"xyz");
  Sty->SetTitleOffset(1.35,"y");
  Sty->SetTitleOffset(1.1,"x");
  Sty->SetTitleFillColor(10);
  Sty->SetTitleTextColor(kBlack);
  Sty->SetTickLength(.03,"xz");
  Sty->SetTickLength(.02,"y");
  Sty->SetNdivisions(210,"x");
  Sty->SetNdivisions(210,"yz");
  Sty->SetEndErrorSize(10);
  gROOT->SetStyle("MyStyle");
  gROOT->ForceStyle();


  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd(1);

  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);

  string histname = "InvMass/6Li/Ex_6Li_da_long";

  TFile *input = new TFile("/home/Li7star/Li7/sort.root", "read");
  TH1I *hist = (TH1I*)input->Get(histname.c_str());

  histname = "Ex_long";
  TFile *input2 = new TFile("/home/Li7star/li7sim/sim_Li6.root", "read");
  TH1I *hist2 = (TH1I*)input2->Get(histname.c_str());


  int Nfits = 1;
  float xlow = 2.13;
  float xhigh = 2.25;

  int nVal = hist->GetNbinsX();
  float x[nVal], y[nVal], a[nVal]; //save the discrete points of
  for (int i=0; i<nVal; i++)
  {
    x[i] = hist->GetBinCenter(i);
    y[i] = hist->GetBinContent(i);
    a[i] = hist2->GetBinContent(i);
    cout << "a[" << i << "] " << a[i] << endl;
  }

  //chi2 *Uncert = new chi2(nVal, hist->FindBin(xlow), hist->FindBin(xhigh), y, a);

  //Fitting functions stored in simulation.
  sle *LinEq = new sle(Nfits);

  //fill "sle" LinEq.Y[k], where y[i] is the solution.
  //k=0 is simulated response
  float sum = 0; 
  for (int i=hist->FindBin(xlow); i<hist->FindBin(xhigh)+1; i++)
  {
    sum += (float)hist->GetBinContent(i) * (float)hist2->GetBinContent(i);
  }
  LinEq->Y[0] = sum;


  //start filling matrix M[Nfits][Nfits] = M[1][1]
  //Get the diagonal elements
  sum = 0;
  for (int i=hist->FindBin(xlow); i<hist->FindBin(xhigh)+1; i++)
  {
    sum += hist2->GetBinContent(i) * hist2->GetBinContent(i);
  }
  LinEq->M[0][0] = (float)sum;
  
  LinEq->solve();

  hist2->Scale(LinEq->Y[0]);


  //calculate uncertainty
  //float chi2val = Uncert->CalcChisqr(LinEq->Y[0]);
  //float sigmaB = Uncert->CalcUncertainty(LinEq->Y[0],LinEq->Y[1],LinEq->Y[2]);

  //cout << "scale sim  " << LinEq->Y[0] << endl;
  //cout << "chi2val    " << chi2val << endl;
  //cout << "scale fep  " << LinEq->Y[1]*1000000 << " ± staterr " << sigmaB*1000000 << endl;
  //cout << "scale comp " << LinEq->Y[2] << endl;
  //cout << "chi2 " << chi2val << endl;

  //fit->SetLineColor(2);
  //fit->SetLineStyle(1);
  //fit->GetXaxis()->SetRangeUser(xlow,xhigh);
  
  //give sqrt(N) as error for each bin in the histogram
  for (int i=0; i<hist->GetNbinsX(); i++)
  {
    hist->SetBinError(i, sqrt(hist->GetBinContent(i)));
  }


  hist->GetXaxis()->SetTitle("E* (MeV)");
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.9);
  hist->SetMarkerStyle(8);
  hist->GetXaxis()->SetRangeUser(1.9,2.5);
  hist->SetStats(kFALSE);


  hist2->SetLineColor(3);
  hist2->SetLineStyle(2);

/*
  //TLegend *legend = new TLegend(0.59,0.8,0.875,0.93);
  //TLegend *legend = new TLegend(0.675,0.77,0.925,0.95);
  TLegend *legend = new TLegend(0.2,0.77,0.40,0.95);
  legend->AddEntry(h_gamma, "Exp, 6Li+#gamma");
  legend->AddEntry(fit, "fitted");
  legend->AddEntry(bkgsmooth, "bkg, Be7+#gamma");
  legend->AddEntry(h_geant_fep, "Geant4 sim, FEP");
  legend->AddEntry(h_geant_comp, "Geant4 sim, Comp");
*/
  hist->Draw("PE");
  //legend->Draw();

  //bkg->Draw("PEhist same");
  //bkgsmooth->Draw("Chist same");

  hist2->Draw("Chist same");
  //h_geant_comp->Draw("Chist same");

  //fit->Draw("Chist same");

/*
  TLatex tt2;
  tt2.SetTextColor(1);
  tt2.SetTextAngle(0.);
  tt2.SetTextSize(tt.GetTextSize()*1.2);
  tt2.DrawLatex(6.1,325,"4125 events");
  ostringstream text;
  text << round(LinEq->Y[1]*1000000) << " #pm " << round(sigmaB*1000000) << endl;
  tt2.DrawLatex(3.8,90,text.str().c_str());
  tt2.DrawLatex(3.8,75,"Gammas emitted");
*/


  string printname = "temp.png";
  c1->Print(printname.c_str(), "png");

  //delete LinEq;
  //delete Uncert;
}
