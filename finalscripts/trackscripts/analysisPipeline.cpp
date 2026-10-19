#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TObject.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <cmath>
void angular_distribution(TTree* cut_tree){
    // Histogram objects 
    TH2D *h1 = new TH2D("h1", "2D Angular distribution", 16, 0, 180, 32, -180, 180);
    TH1D *h2 = new TH1D("h2", "#theta", 100, -180, 180);
    TH1D *h3 = new TH1D("h3", "#phi", 100, -180, 180);
    // Set up variables for vertex coordinates, theta, phi
    // First vertex coords
    std::vector<double>* first_vertex_x = 0;
    std::vector<double>* first_vertex_y = 0;
    std::vector<double>* first_vertex_z = 0;
    // Second vertex coords
    std::vector<double>* second_vertex_x = 0;
    std::vector<double>* second_vertex_y = 0;
    std::vector<double>* second_vertex_z = 0;
    // Zenith angle theta
    double theta = 0;
    // Azimuth angle phi
    double phi = 0;
    cut_tree->SetBranchAddress("first_vertex_x", &first_vertex_x);
    cut_tree->SetBranchAddress("first_vertex_y", &first_vertex_y);
    cut_tree->SetBranchAddress("first_vertex_z", &first_vertex_z);
    cut_tree->SetBranchAddress("second_vertex_x", &second_vertex_x);
    cut_tree->SetBranchAddress("second_vertex_y", &second_vertex_y);
    cut_tree->SetBranchAddress("second_vertex_z", &second_vertex_z);
    // Add new branches for theta and phi
    TBranch *b1 = cut_tree->Branch("theta", &theta);
    TBranch *b2 = cut_tree->Branch("phi", &phi);
    // Event loop
    Long64_t n = cut_tree->GetEntries();
    // Debugging
    // std::cout << n << std::endl;
    for (Long64_t i = 0; i < n; i++){
        cut_tree->GetEntry(i);
        for (Long64_t j = 0; j < first_vertex_x->size(); j++){
            double x = first_vertex_x->at(j) - second_vertex_x->at(j);
            double y = first_vertex_y->at(j) - second_vertex_y->at(j);
            double z = first_vertex_z->at(j) - second_vertex_z->at(j);
            double r = sqrt(x*x + y*y + z*z);
            if (r == 0) continue;
            theta = acos(z/r)*(180/TMath::Pi());
            phi = atan2(y, x)*(180/TMath::Pi());
            // std::cout << phi << std::endl;
            h1->Fill(theta, phi);
            h2->Fill(theta);
            h3->Fill(phi); 
        }
    }
    // Figure window
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
    h1->GetYaxis()->SetTitleSize(0.05);
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
    double thetaerr = h2->GetStdDev();
    std::cout << "Zenith angle mode: " << maxTheta <<  " degrees" << std::endl;
    std::cout << "Mean zenith angle: " << meanTheta << " +/- " << thetaerr << " degrees" << std::endl;
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
    double phierr = h3->GetStdDev();
    std::cout << "Azimuth angle mode: " << maxPhi << " degrees" << std::endl;
    std::cout << "Mean azimuth angle: " << meanPhi << " +/- " << phierr <<" degrees" << std::endl;
    c1->SaveAs("mutrackdistribution2.png");
    TFile *fout = new TFile("trackAnalysis.root", "RECREATE");
    TTree *new_tree = cut_tree->CloneTree(0);
    new_tree->CopyEntries(cut_tree);
    h1->Write();
    h2->Write();
    h3->Write();
    new_tree->Write();
    fout->Close();
    
}

void analysisPipeline(const char* filename){
    // Open file and retrieve tree
    TFile *f = TFile::Open(filename);
    TTree *tree = (TTree*)f->Get("Result_tree");
    // Set up variables for existing branches
    std::vector<double>* energy = nullptr;
    std::vector<int>* om_number = nullptr;
    std::vector<int>* pid = nullptr;
    std::vector<double>* first_vertex_x = nullptr;
    std::vector<double>* second_vertex_x = nullptr;
    std::vector<double>* first_vertex_y = nullptr;
    std::vector<double>* second_vertex_y = nullptr;
    // Set branch addresses
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("first_vertex_x", &first_vertex_x);
    tree->SetBranchAddress("second_vertex_x", &second_vertex_x);
    tree->SetBranchAddress("first_vertex_y", &first_vertex_y);
    tree->SetBranchAddress("second_vertex_y", &second_vertex_y);
    // Make intermediate tree
    gROOT->cd();
    TTree *itree = tree->CloneTree(0);
    // Set up new variable for number of OM hits
    int OM_count = 0;
    TBranch *b = itree->Branch("OM_count", &OM_count);
    // Event loop
    Long64_t m = tree->GetEntries();
    for (Long64_t i = 0; i < m; i++){
        tree->GetEntry(i);
        // Calculate the number of OMs hit 
        OM_count = om_number->size();
        // if (!first_vertex_x || !second_vertex_x || first_vertex_x->size() != second_vertex_x->size()) continue;
        // std::cout << OM_count << std::endl;
        for (int j = 0; j < first_vertex_y->size(); j++){
            for (int k = 0; k < second_vertex_y->size(); k++){
                // if (pid->at(j) == 1 && pid->at(k) == 1){
                //     // std::cout << "1: " << first_vertex_y->at(j) << ": PID: " << pid->at(j) << std::endl;
                //     // std::cout << "2: " << second_vertex_y->at(k) << ": PID: " << pid->at(k) << std::endl;
                // }
            }
            
        }
        // Fill new tree
        itree->Fill();
    }
    // Apply selection cut (E >= 10 MeV, OM_count >= 2, electron_number == 1 for track selection)
    // oops changed cut to 12 MeV 
    // PID: electron is 1 NOT 0
    TTree *cut_tree = itree->CopyTree("energy >= 12 && OM_count >= 2 && electron_number == 2");
    std::cout << cut_tree->GetEntries() << std::endl;
    // Pass cut tree to angular distribution analysis
    angular_distribution(cut_tree);
    
    // Memory clean-up
    delete cut_tree;
    delete itree;
    f->Close();
    
    
}