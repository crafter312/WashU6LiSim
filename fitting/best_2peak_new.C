{
 gROOT->Reset();
  TStyle * Sty = new TStyle("MyStyle","MyStyle");
  Sty->SetOptTitle(0);
  Sty->SetOptStat(0);
  Sty->SetLineWidth(5);
  Sty->SetPalette(55);
  Sty->SetCanvasColor(10);
  Sty->SetCanvasBorderMode(0);
  Sty->SetFrameLineWidth(3);
  Sty->SetFrameFillColor(10);
  Sty->SetPadColor(10);
  Sty->SetPadTickX(1);
  Sty->SetPadTickY(0);
  Sty->SetPadBottomMargin(.2);
  Sty->SetPadTopMargin(.05);
  Sty->SetPadLeftMargin(.13);
  Sty->SetPadRightMargin(.09);
  Sty->SetHistLineWidth(3);
  Sty->SetFuncWidth(3);
  Sty->SetFuncColor(kGreen);
  Sty->SetLineWidth(3);
  Sty->SetLabelSize(0.06,"xyz");
  Sty->SetLabelOffset(0.01,"y");
  Sty->SetLabelOffset(0.02,"x");
  Sty->SetLabelColor(kBlack,"xyz");
  Sty->SetTitleSize(0.07,"xyz");
  Sty->SetTitleOffset(1.1,"y");
  Sty->SetTitleOffset(1.2,"x");
  Sty->SetTitleFillColor(10);
  Sty->SetTitleTextColor(kBlack);
  Sty->SetTickLength(.05,"xz");
  Sty->SetTickLength(.025,"y");
  Sty->SetNdivisions(10,"xyz");
  Sty->SetEndErrorSize(0);
  gROOT->SetStyle("MyStyle");
  gROOT->ForceStyle();

  TCanvas canvas("Ex_ta");
  ifstream file("out.dat");
  double a,b,c,d,e,f,g,h;

  int n = 0;
  double x[500];
  double y[500];
  double sx[500];
  double sy[500];
  double ytot[500];
  double yb[500];
  double yy1[500];
  double yy2[500];

  //sim of new state
  double yy3[500];
  double yy3a[500];
  double yy4[500];
  double yy5[500];
  double yy6[500];

  for (;;)
  {
    file >> a >> b >> c>> d >> e >> f >> g;
    if (file.eof())break;
    if (file.bad())break;

    x[n] = a;
    y[n] = b;
    sy[n] = c;
    sx[n] = 0.;
    ytot[n] = d;
    yb[n] = e;

    yy1[n] = f;
    yy2[n] = g;
    //yy6[n] = h;
    n++;
  }

  //need to get simulation of 
  TFile * fsim_new_pureP = new TFile ("out_sim_10205_130.root");
  TH1I* sim_new_pureP = (TH1I*) fsim_new_pureP->Get("Ex");
  sim_new_pureP->Rebin(2);
  sim_new_pureP->Smooth(2);

  TFile * fsim_new_pureP2 = new TFile ("../rootout/sim_p6He_10205_130.root");
  TH1I* sim_new_pureP2 = (TH1I*) fsim_new_pureP2->Get("Ex");
  sim_new_pureP2->Rebin(2);
  sim_new_pureP2->Smooth(2);

  TFile * fsim_new_mixed = new TFile ("../rootout/sim_line_10200_130_longline.root");
  TH1I* sim_new_mixed = (TH1I*) fsim_new_mixed->Get("Ex");
  sim_new_mixed->Rebin(2);
  sim_new_mixed->Smooth(2);

  TFile * fsim_new_mixedhard = new TFile ("../rootout/sim_line_10040_130_longline.root");
  TH1I* sim_new_mixedhard = (TH1I*) fsim_new_mixedhard->Get("Ex");
  sim_new_mixedhard->Rebin(2);
  sim_new_mixedhard->Smooth(2);


  TFile * file_eff = new TFile ("../rootout/sim_uniform.root");
  TH1I* hist_eff = (TH1I*) file_eff->Get("Ex");

  hist_eff->Rebin(2);

  double eff[400];
  double eff_norm[400];
  double x_eff[400];

  double max = 0;
  //find max value in response to uniform energy
  for (int i=0;i<400;i++)
  {
    if (hist_eff->GetBinContent(i) > max) max = hist_eff->GetBinContent(i);
  }
  //normalize efficiency to max value
  for (int i=0;i<120;i++)
  {
    if (hist_eff->GetBinContent(i) > 0)
    {
      eff[i] = hist_eff->GetBinContent(i)/max*700;
      eff_norm[i] = hist_eff->GetBinContent(i)/max;
      x_eff[i] = hist_eff->GetBinCenter(i);
      cout << "i" << i << " x " << x_eff[i] << "  eff " << eff[i] << endl;
    }
    else
    {
      eff[i] = 0;
      eff_norm[i] = 0;
    }
  }

  hist_eff->Smooth(3);


  xlow = 10.0;
  xhigh = 13.2;
  int low = sim_new_pureP->FindBin(xlow);
  int high = sim_new_pureP->FindBin(xhigh);

  cout << "n " << n << " vs high-low " << high << " - " << low << endl;

  for (int j=low;j<high;j++)
  {
    int i = j - low; 


    //(from para2.dat)*(fit.cpp ln372) * (15 mil in IAS sim vs 2 mil in new sim) * (BR and C2S factor) * (scale)
    //yy3[i] = sim_new_pureP->GetBinContent(j) * 1659.03 * 1e-6 * (15/2) * (1/0.07)*0.02;

     //(from para2.dat)*(fit.cpp ln372) * (15 mil in IAS sim vs 2 mil in new sim) * (BR and C2S factor of populating IAS) * (specroscopic factor of this state)
    double scale = 25574.7 * 1e-6 * (1./2.) * (1/0.07)*0.02;
    yy3a[i] = sim_new_pureP2->GetBinContent(j)/eff_norm[j] * scale;
    //cout << "yy3[" << i << "] = " << yy3[i] << endl;
    //cout << "yy3a[" << i << "] = " << yy3a[i] << endl;


     //(from para2.dat)*(fit.cpp ln372) * (1 mil in IAS sim vs 2 mil in new sim) * (BR and C2S factor of populating IAS) * (branching ratio p/n * specroscopic factor of this state) * (fresco ratio of xsec(E) )
    scale = 25574.7 * 1e-6 * (1./20.) * (1/0.07) * 0.583 * 0.9 * ((-12.32*sim_new_mixed->GetBinCenter(j) + 214.41)/209.75) /eff_norm[j];
    //cout << "E = " << sim_new_mixed->GetBinCenter(j) << "  ratio " << (-12.32*sim_new_mixed->GetBinCenter(j) + 214.41)/209.75 << "  scale " << scale << endl;
    yy4[i] = sim_new_mixed->GetBinContent(j) * scale;


    scale = 25574.7 * 1e-6 * (1./20.) * (1/0.07) * 0.45 * 0.9 * ((-12.32*sim_new_mixedhard->GetBinCenter(j) + 214.41)/209.75) /eff_norm[j];
    yy5[i] = sim_new_mixedhard->GetBinContent(j) * scale;

  }

/*
  low = sim_new_pureP->FindBin(xlow);
  high = sim_new_pureP->FindBin(xhigh);

  cout << "n " << n << " vs high-low " << high << " - " << low << endl;

  for (int j=low;j<high;j++)
  {
    int i = j - low; 
    //(from para2.dat)*(fit.cpp ln372) * (15 mil in IAS sim vs 2 mil in new sim) * (BR and C2S factor) * (scale)
    yy3[i] = sim_new_pureP->GetBinContent(j) * 1659.03 * 1e-6 * (15/2) * (1/0.07)*0.02;
    yy3a[i] = sim_new_pureP2->GetBinContent(j)/eff_norm[j]  * 1659.03 * 1e-6 * (15/2) * (1/0.07)*0.02;
    cout << "yy3[" << i << "] = " << yy3[i] << endl;
    cout << "yy3a[" << i << "] = " << yy3a[i] << endl;
  }

*/




  TH2S frame("frame","",10,9.98,13.0,10,0,750);
  //TH2S frame("frame","",10,9.98,12.0,10,0,180);
  frame.GetXaxis()->SetTitle("E* (MeV)");
  frame.GetYaxis()->SetTitle("Counts / 40 keV");
  frame.GetYaxis()->SetTitleSize(0.06);
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.Draw();

  TGraphErrors gg(n,x,y,sx,sy);
  gg.SetMarkerStyle(20);
  gg.SetMarkerColor(1);
  gg.SetLineColor(1);
  gg.SetMarkerSize(1.5);
  gg.Draw("PE");

  TGraph gtot(n,x,ytot);
  gtot.SetLineColor(2);


  //state2
  TGraph gb(n,x,yb);
  gb.SetLineColor(4);
  gb.SetLineStyle(9);
  gb.Draw("C");

  //state1
  TGraph g1(n,x,yy1);
  g1.SetLineColor(8);
  g1.SetLineStyle(2);
  g1.SetLineWidth(4);
  g1.Draw("C");

  //state2
  TGraph g2(n,x,yy2);
  g2.SetLineColor(8);
  g2.SetLineStyle(2);
  g2.SetLineWidth(4);
  g2.Draw("C");

  //predicted did in the old way of doing it
  //TGraph g3(n,x,yy3);
  //g3.SetLineColor(6);
  //g3.SetLineStyle(2);
  //g3.SetLineWidth(4);
  //g3.Draw("C");

  //predicted_a
  TGraph g3a(n,x,yy3a);
  g3a.SetLineColor(6); //magenta
  g3a.SetLineStyle(3);
  g3a.SetLineWidth(5);
  g3a.Draw("C");

  //predicted but lineshape comes from sim
  TGraph g4(n,x,yy4);
  g4.SetLineColor(7); //cyan
  g4.SetLineStyle(3);
  g4.SetLineWidth(5);
  g4.Draw("C");

  //predicted but lineshape comes from sim
  TGraph g5(n,x,yy5);
  g5.SetLineColor(kOrange-3);
  g5.SetLineStyle(3);
  g5.SetLineWidth(5);
  g5.Draw("C");

  //Fitted sim state to show scale
  TGraph g6(n,x,yy6);
  g6.SetLineColor(7);
  g6.SetLineStyle(2);
  g6.SetLineWidth(4);
  //g6.Draw("C");

  TGraph geff(100,x_eff,eff);
  geff.SetLineColor(14);
  geff.SetLineStyle(9);
  geff.Draw("L");

  //move the toal plot to the end
  gtot.Draw("L");


  TGaxis *axis = new TGaxis(13, 0, 13, 700,
                           0,1,5,"+L");
  axis->SetLineColor(14);
  axis->SetLabelColor(14);

  axis->SetTitleColor(14);
  axis->SetTitleSize(0.06);
  axis->SetTitle("relative efficiency");
  axis->SetTitleOffset(0.65);
  axis->CenterTitle();
  axis->Draw();




  TLatex tt;
  tt.SetTextSize(.07);
  tt.DrawLatex(11.75, 550,"^{7}Li#rightarrow#it{p}+^{6}He");
  TArrow arrow;
  arrow.SetFillColor(1);
  arrow.SetAngle(35);

  //arrow.DrawArrow(4.652,260,4.652,230,.02);
  //arrow.DrawArrow(7.454,60,7.454,30,.02);
  //arrow.DrawArrow(3.05,650,3.05,550,.02);
  //arrow.DrawArrow(4.35,500,4.35,400,.02);

  //tt.DrawLatex(4.7,260,"#frac{7}{2}^{-}");
  //tt.DrawLatex(7.55,60,"#frac{5}{2}^{-}");
  //tt.DrawLatex(2.,512,"1^{-}");
  //tt.DrawLatex(1.42,512,"2^{-}");

  TLegend *legend = new TLegend(0.64,0.58,0.88,0.85);
  legend->AddEntry(&gg, "Data","ep");
  legend->AddEntry(&gtot, "Fit");
  legend->AddEntry(&g1, "Sim Peaks");
  legend->AddEntry(&gb, "Linear Bkg");
  legend->AddEntry(&g3a, "Near Threshold Res");
  //legend->Draw();



 
  string printname = "Ex_7Li_p6He_temp.png";
  canvas.Print(printname.c_str(), "png");
  printname = "Ex_7Li_p6He_temp.eps";
  canvas.Print(printname.c_str(), "eps");
  printname = "Ex_7Li_p6He_temp.pdf";
  canvas.Print(printname.c_str(), "pdf");
}
