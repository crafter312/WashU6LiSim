#include <TFile.h>
#include <TTree.h>

void rotationalias() {
    TFile* ifile = TFile::Open("sim_dalpha_3+_56MeV_90mm_diamond100keV_loss0-005step_0recoilquench.root");
    TTree* t = (TTree*)ifile->Get("t");

		// 1. Force ROOT to drop any previous memory associations
    t->SetBranchStatus("*", 1); 
    t->ResetBranchAddresses();

    // Print the first few entries of the ACTUAL branch to the terminal
    std::cout << "Checking actual branch data:" << std::endl;
    t->Scan("targEloss", "isFragDet", "", 5); 

    t->SetAlias("y", "targEloss");
    
    std::cout << "Checking alias data:" << std::endl;
    t->Scan("y", "isFragDet", "", 5);
    
    t->Draw("y");
}
