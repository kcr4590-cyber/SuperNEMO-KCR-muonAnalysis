#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <queue>
#include <map>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TObject.h>
#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TSystemDirectory.h>
#include <TSystem.h>
#include <TList.h>

void trackAnalysis(){
    TString parent = "/sps/nemo/scratch/Edinburgh_student/Katherine/TRACK2/"; // path to folder containing individual run result folders
    TString prefix = "run-"; // prefix of every directory name
    TString filename = "trackAnalysis.root"; // common file name
    TString treeName = "Result_tree";
    TString hist_name = "h1";
    TString hist_name2 = "h2";
    TString hist_name3 = "h3";


    // Make a TChain to link all the files together
    TChain* chain = new TChain(treeName);

    // Look inside parent directory
    TSystemDirectory sys_dir(parent, parent);
    TList* filelist = sys_dir.GetListOfFiles();
    if (!filelist){
        std::cerr << "Can't read directory " << parent << std::endl;
        return;
    }

    // Loop through file list to find matching directories
    TIter next(filelist);
    TSystemFile* file;
    while ((file = (TSystemFile*)next())){
        TString name = file->GetName();
        // Check if item is directory and matches prefix
        if (file->IsDirectory() && name.BeginsWith(prefix) && name != "." && name != ".."){
            TString path = parent + name + "/" + filename;

            if (!gSystem->AccessPathName(path)){
                chain->Add(path); // Add the file to the TChain
                std::cout << "Added: " << path << std::endl;
            }
            else{
                std::cout << "File missing in folder: " << path << std::endl;
            }
        }
    }
    // Get data from trees in the chain
    double theta;
    double phi;
    chain->SetBranchAddress("theta", &theta);
    chain->SetBranchAddress("phi", &phi);
    // Output file
    TFile *f = new TFile("trackhistograms.root", "RECREATE");
    TH2D *h1 = new TH2D("h1", "2D Angular distribution", 16, 0, 180, 32, -180, 180);
    h1->SetCanExtend(TH2::kAllAxes);
    TH1D *h2 = new TH1D("h2", "#theta", 100, -180, 180);
    h2->SetCanExtend(TH1::kAllAxes);
    TH1D *h3 = new TH1D("h3", "#phi", 100, -180, 180);
    h3->SetCanExtend(TH1::kAllAxes);
    Long64_t n = chain->GetEntries();
    std::cout << n << std::endl;
    int currentTree = -1;
    for (Long64_t i = 0; i < n; i++){
        Long64_t entry = chain->LoadTree(i);
        if (entry < 0) break;
        if (chain->GetTreeNumber() != currentTree){
            currentTree = chain->GetTreeNumber();
            chain->SetBranchAddress("theta", &theta);
            chain->SetBranchAddress("phi", &phi);
        }

        chain->GetEntry(i);
        h1->Fill(theta, phi);
        h2->Fill(theta);
        h3->Fill(phi);
    }

    TCanvas *c1 = new TCanvas("c1", "Angular Distributions", 1600, 1200);
    c1->Divide(2,2);
    c1->cd(1);
    // Draw the histogram
    h1->SetStats(0);
    h1->Draw("COLZ");
    gStyle->SetPalette(kCool);
    h1->GetXaxis()->SetTitle("Elevation angle #theta (degrees)");
    h1->GetYaxis()->SetTitle("Azimuthal angle #phi (degrees)");
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleOffset(0.95);
    h1->GetXaxis()->SetLabelSize(0.04);
    h1->GetYaxis()->SetLabelSize(0.04);
    c1->cd(2);
    h2->SetStats(0);
    h2->Draw();
    h2->GetXaxis()->SetTitle("Elevation angle #theta (degrees)");
    h2->GetYaxis()->SetTitle("Entries");
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->GetXaxis()->SetLabelSize(0.04);
    h2->GetYaxis()->SetLabelSize(0.04);
    double meanTheta = h2->GetMean();
    int binmax = h2->GetMaximumBin();
    double maxTheta = h2->GetXaxis()->GetBinCenter(binmax);
    double maxThetaerr = h2->GetBinError(binmax);
    double thetaerr = h2->GetStdDev();
    std::cout << "Zenith angle mode: " << maxTheta << " +/- " << maxThetaerr <<  " degrees" << std::endl;
    std::cout << "Mean zenith angle: " << meanTheta << " +/- " << thetaerr << " degrees" << std::endl;
    TF1 *fit = new TF1("fit", "gaus(0) + gaus(3)", 0, 180);
    fit->SetParameter(0, 30);
    fit->SetParameter(1, 44);
    fit->SetParameter(2, 8);
    fit->SetParameter(3, 20);
    fit->SetParameter(4, 140);
    fit->SetParameter(5, 8);
    fit->SetParLimits(1, 25, 90);
    fit->SetParLimits(4, 120, 180);
    double max_sigma = 8;
    fit->SetParLimits(2, 0.1, max_sigma);
    fit->SetParLimits(5, 0.1, max_sigma);
    h2->Fit(fit, "RB");
    fit->SetLineColor(kBlack);
    fit->Draw("SAME");
    c1->cd(3);
    h3->SetStats(0);
    h3->Draw();
    h3->GetXaxis()->SetTitle("Azimuth angle #phi (degrees)");
    h3->GetYaxis()->SetTitle("Entries");
    h3->GetXaxis()->SetTitleSize(0.05);
    h3->GetYaxis()->SetTitleSize(0.05);
    h3->GetXaxis()->SetLabelSize(0.04);
    h3->GetYaxis()->SetLabelSize(0.04);
    double meanPhi = h3->GetMean();
    int binmax2 = h3->GetMaximumBin();
    double maxPhi = h3->GetXaxis()->GetBinCenter(binmax2);
    double maxPhierr = h3->GetBinError(binmax2);
    double phierr = h3->GetStdDev();
    std::cout << "Azimuth angle mode: " << maxPhi << " +/- " << maxPhierr << " degrees" << std::endl;
    std::cout << "Mean azimuth angle: " << meanPhi << " +/- " << phierr <<" degrees" << std::endl;
    c1->SaveAs("muangulardist2.png");
    h1->Write();
    h2->Write();
    h3->Write();
    f->Close();

}
