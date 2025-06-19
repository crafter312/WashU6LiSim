{
	gROOT->Reset();

	// Set default style attributes
	TStyle * Sty = new TStyle("MyStyle","MyStyle");
	Sty->SetOptTitle(0);
	Sty->SetOptStat(0);
	Sty->SetLineWidth(4);
	gStyle->SetPalette(kBird);
	Sty->SetCanvasColor(10);
	Sty->SetCanvasBorderMode(0);
	Sty->SetFrameLineWidth(0);
	Sty->SetFrameFillColor(10);
	Sty->SetPadColor(10);
	Sty->SetPadTickX(1);
	Sty->SetPadTickY(1);
	Sty->SetPadBottomMargin(.15);
	Sty->SetPadTopMargin(.03);
	Sty->SetPadLeftMargin(.15);
	Sty->SetPadRightMargin(.15);
	Sty->SetHistLineWidth(3);
	Sty->SetFuncWidth(3);
	Sty->SetFuncColor(kGreen);
	Sty->SetLineWidth(3);
	Sty->SetLabelSize(0.05,"xyz");
	Sty->SetLabelOffset(0.015,"y");
	Sty->SetLabelOffset(0.02,"x");
	Sty->SetLabelColor(kBlack,"xyz");
	Sty->SetTitleSize(0.06,"y");
	Sty->SetTitleSize(0.07,"x");
	Sty->SetTitleOffset(0.95,"y");
	Sty->SetTitleOffset(1.2,"x");
	Sty->SetTitleFillColor(10);
	Sty->SetTitleTextColor(kBlack);
	Sty->SetTickLength(.05,"xz");
	Sty->SetTickLength(.025,"y");
	Sty->SetNdivisions(10,"y");
	Sty->SetNdivisions(10,"x");
	Sty->SetEndErrorSize(0);
	Sty->SetTextFont(42);
	gROOT->SetStyle("MyStyle");
	gROOT->ForceStyle();

	// Commonly modified style attributes
	Sty->SetNdivisions(10,"xy");  // # tick mark divisions on the x and y axes
	Sty->SetTitleOffset(0.7,"y"); // y axis title offset
	Sty->SetTitleSize(0.06,"y");  // y axis title size
	Sty->SetTitleOffset(1.1,"x"); // x axis title offset
	Sty->SetTitleSize(0.06,"x");  // x axis title size

	Sty->SetPadBottomMargin(.15);
	Sty->SetPadTopMargin(.08);
	Sty->SetPadLeftMargin(.12);
	Sty->SetPadRightMargin(.12);

	/******** HISTOGRAM SETUP ********/

	TCanvas c1("c1","",1202,576);

	TFile* ifile = TFile::Open("../rootfiles/sim_alphapn_56MeV_90mm.root");
	TTree* t = (TTree*)ifile->Get("t");
	TH2I kin("kin","",120,0,60,40,0,20);
	t->Draw("ENeut:thetaNeut>>kin","","goff");

	// Set graphical properties
	kin.GetXaxis()->SetTitle("Neutron polar angle (degrees)");
	kin.GetXaxis()->CenterTitle();
	kin.GetYaxis()->SetTitle("Neutron energy (MeV)");
	kin.GetYaxis()->CenterTitle();
	kin.Draw("colz");

	// Save figure
	c1.Print("6Lisim_EvsThetaNeut_56MeV90mm.png", "png");
	c1.Print("6Lisim_EvsThetaNeut_56MeV90mm.eps", "eps");
}



