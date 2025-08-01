void plot_sim(){
  
  TFile *input = new TFile("sim.root", "read");
  
  bool save = true;
  
  bool plotopt[4] = {1,0,0,0};
 
 
  ////////////////////////////////////////////////////////////////////////////
  // Canvas#1
  ////////////////////////////////////////////////////////////////////////////
  if (plotopt[0]){ 
    TCanvas *c1 = new TCanvas("c1", "vel"); //, 1600, 1200);
    c1->SetCanvasSize(1600,1200);
    c1->SetWindowSize(800,600);
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetTickx();
    gPad->SetTicky();

    TH1F *hist_vel_P = (TH1F*)input->Get("vel_P");
    hist_vel_P->SetStats(0);
    //hist_vel_P->SetTitle("Velocity distribution of Li7* in Lab");
    hist_vel_P->GetXaxis()->SetTitle("Velocity (cm/ns)");
    hist_vel_P->GetYaxis()->SetRangeUser(0,50000);
    hist_vel_P->GetXaxis()->CenterTitle();
    hist_vel_P->GetXaxis()->SetTitleSize(0.05);
    hist_vel_P->GetXaxis()->SetLabelSize(0.05);
    hist_vel_P->Draw();


    TH1F *hist_vel_S = (TH1F*)input->Get("vel_S");
    hist_vel_S->SetStats(0);
    hist_vel_S->SetLineColor(kRed);
    hist_vel_S->Draw("same");
    
    TH1F *hist_vel_T = (TH1F*)input->Get("vel_T");
    hist_vel_T->SetStats(0);
    hist_vel_T->SetLineColor(kGreen);
    hist_vel_T->Draw("same");

    //if (save) c1->SaveAs("plotting/vel.jpg");
    
    ////////////////////////////////////////////////////////////////////////////
    c1->cd(2);
    gPad->SetTickx();
    gPad->SetTicky();

    TH1F *hist_theta_P = (TH1F*)input->Get("theta_P");
    hist_theta_P->SetTitle("Theta districtuion of Li7* in Lab");
    hist_theta_P->SetStats(0);
    hist_theta_P->GetXaxis()->SetTitle("Theta (deg)");
    hist_theta_P->GetXaxis()->CenterTitle();
    hist_theta_P->GetXaxis()->SetTitleSize(0.05);
    hist_theta_P->GetXaxis()->SetLabelSize(0.05);
    hist_theta_P->Draw();


    TH1F *hist_theta_S = (TH1F*)input->Get("theta_S");
    hist_theta_S->SetStats(0);
    hist_theta_S->SetLineColor(kRed);
    hist_theta_S->Draw("same");
    
    TH1F *hist_theta_T = (TH1F*)input->Get("theta_T");
    hist_theta_T->SetStats(0);
    hist_theta_T->SetLineColor(kGreen);
    hist_theta_T->Draw("same");

    //if (save) c2->SaveAs("plotting/theta.jpg");

    ////////////////////////////////////////////////////////////////////////////
    c1->cd(3);
    gPad->SetTickx();
    gPad->SetTicky();

    TH1F *hist_phi_P = (TH1F*)input->Get("phi_P");
    hist_phi_P->SetTitle("Phi districtuion of Li7* in Lab");
    hist_phi_P->GetYaxis()->SetRangeUser(0,5900);
    hist_phi_P->SetStats(0);
    hist_phi_P->GetXaxis()->SetTitle("Phi (deg)");
    hist_phi_P->GetXaxis()->CenterTitle();
    hist_phi_P->GetXaxis()->SetTitleSize(0.05);
    hist_phi_P->GetXaxis()->SetLabelSize(0.05);
    hist_phi_P->Draw();


    TH1F *hist_phi_S = (TH1F*)input->Get("phi_S");
    hist_phi_S->SetStats(0);
    hist_phi_S->SetLineColor(kRed);
    hist_phi_S->Draw("same");
    
    TH1F *hist_phi_T = (TH1F*)input->Get("phi_T");
    hist_phi_T->SetStats(0);
    hist_phi_T->SetLineColor(kGreen);
    hist_phi_T->Draw("same");

    //if (save) c3->SaveAs("plotting/phi.jpg");

    ////////////////////////////////////////////////////////////////////////////
    c1->cd(4);
    gPad->SetTickx();
    gPad->SetTicky();

    TH1F *hist_theta_neut_P = (TH1F*)input->Get("theta_neut_P");
    hist_theta_neut_P->SetTitle("Theta districtuion of neutron from d(He6,Li7)Neut in Lab");
    hist_theta_neut_P->SetStats(0);
    hist_theta_neut_P->GetXaxis()->SetTitle("Theta (deg)");
    hist_theta_neut_P->GetXaxis()->CenterTitle();
    hist_theta_neut_P->GetXaxis()->SetTitleSize(0.05);
    hist_theta_neut_P->GetXaxis()->SetLabelSize(0.05);
    hist_theta_neut_P->Draw();


    TH1F *hist_theta_neut_S = (TH1F*)input->Get("theta_neut_S");
    hist_theta_neut_S->SetStats(0);
    hist_theta_neut_S->SetLineColor(kRed);
    hist_theta_neut_S->Draw("same");
    
    TH1F *hist_theta_neut_T = (TH1F*)input->Get("theta_neut_T");
    hist_theta_neut_T->SetStats(0);
    hist_theta_neut_T->SetLineColor(kGreen);
    hist_theta_neut_T->Draw("same");

    //if (save) c4->SaveAs("plotting/theta_neut.jpg");
    if (save) c1->Print("plotting/mutli.png", "png");
    if (save) c1->Print("plotting/mutli.svg", "svg");
  }
  ////////////////////////////////////////////////////////////////////////////
  // Canvas#2
  ////////////////////////////////////////////////////////////////////////////
  if (plotopt[1]){
    TCanvas *c2 = new TCanvas("c2", "kinematic_circle"); //, 1600, 1200);
    if (save) c2->SetWindowSize(1600,1200);
    c2->cd();
    c2->SetLogz();
    gPad->SetTickx();
    gPad->SetTicky();

    TH1F *kinematic_circle = (TH1F*)input->Get("kinematic_circle");

    kinematic_circle->SetStats(0);
    kinematic_circle->SetTitle("Kinematic circle for Li7");
    kinematic_circle->GetXaxis()->SetTitle("X-Velocity (cm/ns)");
    kinematic_circle->GetXaxis()->CenterTitle();
    kinematic_circle->GetXaxis()->SetTitleSize(0.04);
    kinematic_circle->GetXaxis()->SetLabelSize(0.04);
    kinematic_circle->GetYaxis()->SetTitle("Z-Velocity (cm/ns)");
    kinematic_circle->GetYaxis()->CenterTitle();
    kinematic_circle->GetYaxis()->SetTitleSize(0.04);
    kinematic_circle->GetYaxis()->SetLabelSize(0.04);
    

    kinematic_circle->Draw("colz");
    
    if (save) c2->Print("plotting/kin_circle.png", "png");
    if (save) c2->Print("plotting/kin_circle.svg", "svg");
  }
  ////////////////////////////////////////////////////////////////////////////
  // Canvas#3
  ////////////////////////////////////////////////////////////////////////////
  if (plotopt[2]){
    TCanvas *c3 = new TCanvas("c3", "Gobi_hitpattern"); //, 1600, 1200);
    if (save) c3->SetWindowSize(1600,600);
    c3->Divide(2,1);
    c3->cd(1);
    gPad->SetTickx();
    gPad->SetTicky();
    
    TH1F *protonXY_S = (TH1F*)input->Get("protonXY_S");

    protonXY_S->SetStats(0);
    protonXY_S->SetTitle("Proton hitpattern");
    protonXY_S->GetXaxis()->SetTitle("X (mm)");
    protonXY_S->GetXaxis()->CenterTitle();
    protonXY_S->GetXaxis()->SetTitleSize(0.04);
    protonXY_S->GetXaxis()->SetLabelSize(0.04);
    protonXY_S->GetYaxis()->SetTitle("Y (mm)");
    protonXY_S->GetYaxis()->CenterTitle();
    protonXY_S->GetYaxis()->SetTitleSize(0.04);
    protonXY_S->GetYaxis()->SetLabelSize(0.04);
    
    protonXY_S->Draw("colz");
    //////////////////////////////////////////////////////////////////////////
    c3->cd(2);

    gPad->SetTickx();
    gPad->SetTicky();

    TH1F *coreXY_S = (TH1F*)input->Get("coreXY_S");

    coreXY_S->SetStats(0);
    coreXY_S->SetTitle("He6 hitpattern");
    coreXY_S->GetXaxis()->SetTitle("X (mm)");
    coreXY_S->GetXaxis()->CenterTitle();
    coreXY_S->GetXaxis()->SetTitleSize(0.04);
    coreXY_S->GetXaxis()->SetLabelSize(0.04);
    coreXY_S->GetYaxis()->SetTitle("Y (mm)");
    coreXY_S->GetYaxis()->CenterTitle();
    coreXY_S->GetYaxis()->SetTitleSize(0.04);
    coreXY_S->GetYaxis()->SetLabelSize(0.04);
    
    coreXY_S->Draw("colz");
    
    if (save) c3->Print("plotting/Gobi_hitpattern.png", "png");
    if (save) c3->Print("plotting/Gobi_hitpattern.svg", "svg");

  }
  ////////////////////////////////////////////////////////////////////////////
  // Canvas#4
  ////////////////////////////////////////////////////////////////////////////
  if (plotopt[3]){
    TCanvas *c4 = new TCanvas("c4", ""); //, 1600, 1200);
    if (save) c4->SetWindowSize(1600,1200);
    c4->cd();
    //c4->SetRightMargin(0.09);
    c4->SetLeftMargin(0.12);
    c4->SetBottomMargin(0.12);
    gPad->SetTickx();
    gPad->SetTicky();
    
    TH1F *hist_neut_E_P = (TH1F*)input->Get("neut_E_P");
    hist_neut_E_P->SetStats(0);
    hist_neut_E_P->SetTitle("Energy distribution of neutron -- d(He6,Li7)n");
    hist_neut_E_P->GetXaxis()->SetTitle("Energy (MeV)");
    hist_neut_E_P->GetYaxis()->SetTitle("counts");
    hist_neut_E_P->GetYaxis()->SetRangeUser(0,6000);
    hist_neut_E_P->GetXaxis()->CenterTitle();
    hist_neut_E_P->GetXaxis()->SetTitleSize(0.05);
    hist_neut_E_P->GetYaxis()->SetTitleSize(0.05);
    hist_neut_E_P->GetXaxis()->SetLabelSize(0.05);
    hist_neut_E_P->Draw();


    TH1F *hist_neut_E_S = (TH1F*)input->Get("neut_E_S");
    hist_neut_E_S->SetStats(0);
    hist_neut_E_S->SetLineColor(kRed);
    hist_neut_E_S->Draw("same");
    
    TH1F *hist_neut_E_T = (TH1F*)input->Get("neut_E_T");
    hist_neut_E_T->SetStats(0);
    hist_neut_E_T->SetLineColor(kGreen);
    hist_neut_E_T->Draw("same");
    
    if (save) c4->Print("plotting/neut_energy.png", "png");
    if (save) c4->Print("plotting/neut_energy.svg", "svg");
  }
}
