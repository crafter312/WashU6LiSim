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
  Sty->SetPadTickY(1);
  Sty->SetPadBottomMargin(.2);
  Sty->SetPadTopMargin(.05);
  Sty->SetPadLeftMargin(.15);
  Sty->SetPadRightMargin(.03);
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
  double a,b,c,d,e,f,g;

  int n = 0;
  double x[500];
  double y[500];
  double sx[500];
  double sy[500];
  double ytot[500];
  double yb[500];
  double yy1[500];
  double yy2[500];
  //double yy3[500];
  //double yy4[500];
  //double yy5[500];


  for (;;)
  {
    file >> a >> b >> c>> d >> e >> f;
    if (file.eof())break;
    if (file.bad())break;

    x[n] = a;
    y[n] = b;
    sy[n] = c;
    sx[n] = 0.;
    ytot[n] = d;
    yb[n] = e;

    yy1[n] = f;
    n++;
  }


  TH2S frame("frame","",10,2.5,6,10,0,1100);
  frame.GetXaxis()->SetTitle("E* (MeV)");
  frame.GetYaxis()->SetTitle("Counts / 40 keV");
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
  gtot.Draw("L");

  TGraph gb(n,x,yb);
  gb.SetLineColor(4);
  gb.SetLineStyle(9);
  gb.Draw("C");


  TGraph g1(n,x,yy1);
  g1.SetLineColor(3);
  g1.SetLineStyle(2);
  g1.Draw("C");



  TLatex tt;
  tt.SetTextSize(.07);
  tt.DrawLatex(6.3, 190,"^{7}Li#rightarrow#it{t}+#alpha");
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
 
  string printname = "Ex_7Li_ta_1peak.png";
  canvas.Print(printname.c_str(), "png");
}
