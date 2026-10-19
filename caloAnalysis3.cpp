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
#include <TRatioPlot.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>


void sortVectors(std::vector<double>& vec, std::vector<double>& tdc){
    // Sorts all vectors by increasing tdc
    // Make a vector of indices and sort them by tdc
    if (vec.empty()|| tdc.empty()) return;
    std::vector<size_t> indx(vec.size());
    std::iota(indx.begin(), indx.end(), 0);
    std::sort(indx.begin(), indx.end(), [&](size_t a, size_t b){
        return tdc.at(a) < tdc.at(b);
    });
    // Containers for sorted vectors
    std::vector<double> sorted_vec;
    std::vector<double> sorted_tdc;
    sorted_vec.reserve(vec.size());
    sorted_tdc.reserve(vec.size());
    // Push back the values of the vector in the order of the index vector
    for (size_t idx : indx){
        sorted_vec.push_back(vec.at(idx));
        sorted_tdc.push_back(tdc.at(idx));
    }
    // Overwrite vectors so stored order remains
    vec = std::move(sorted_vec);
    tdc = std::move(sorted_tdc);
    // Sanity check
    // for (size_t i = 0; i < vec.size(); i++){
    //     std::cout << i << (vec)[i] << std::endl;
    //     std::cout << i << (tdc)[i] << std::endl;
    // }
}

void sortintVec(std::vector<int>& vec, std::vector<double>& tdc){
    if (vec.empty() || tdc.empty()) return;
    std::vector<size_t> indx(vec.size());
    std::iota(indx.begin(), indx.end(), 0);
    std::sort(indx.begin(), indx.end(), [&](size_t a, size_t b){
	return tdc.at(a) < tdc.at(b);
    });
    std::vector<int> sorted_vec;
    std::vector<double> sorted_tdc;
    sorted_vec.reserve(vec.size());
    sorted_tdc.reserve(vec.size());
    for (size_t idx: indx){
	sorted_vec.push_back(vec.at(idx));
	sorted_tdc.push_back(tdc.at(idx));
    }
    vec = std::move(sorted_vec);
    tdc = std::move(sorted_tdc);
}    
void sortallVectors(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
std::vector<double>& energy, std::vector<double>& tdc, std::vector<int>& om){
    // Just sorts vectors. idk what else to say about it
    sortVectors(x, tdc);
    sortVectors(y, tdc);
    sortVectors(z, tdc);
    sortVectors(energy, tdc);
    sortintVec(om, tdc);
    std::cout << "successfully sorted vectors" << std::endl;
}

// void checkDirection(const std::vector<double>* vec, int& p2n, int& n2p){
//     p2n = 0;
//     n2p = 0;
//     if (!vec || vec->size() < 2) return;
//     for (size_t i = 0; i < vec->size() - 1; i++){
//         double current = (*vec)[i];
//         double next = (*vec)[i+1];
//         if (current > 0 && next < 0){
//             p2n++;
//         }
//         else if (current < 0 && next > 0){
//             n2p++;
//         }

//     }
// }

// void checkCrossing(const std::vector<double>* x, const std::vector<double>* y, const std::vector<double>* z,
// int& x_p2n, int& x_n2p, int& y_p2n, int& y_n2p, int& z_p2n, int& z_n2p){
//     checkDirection(x, x_p2n, x_n2p);
//     checkDirection(y, y_p2n, y_n2p);
//     checkDirection(z, z_p2n, z_n2p);
// }
void cutClusters(const std::vector<double>* e, const std::vector<double>* t, const std::vector<double>* x,
const std::vector<double>* y, const std::vector<double>* z, const std::vector<int>* om, std::vector<double>& eout,
std::vector<double>& tout, std::vector<double>& xout, std::vector<double>& yout, std::vector<double>& zout, 
std::vector<int>& omout, double max_dt = 20.0){
    eout.clear();
    tout.clear();
    xout.clear();
    yout.clear();
    zout.clear();
    omout.clear();
    std::cout << "made it to cut func" << std::endl;
    std::cout << e->size() << std::endl;
    if (e->empty()|| t->empty() || x->empty() || y->empty() || z->empty()) return;
    // Sort by highest energy
    std::vector<int> indices(e->size());
    for (size_t i = 0; i < indices.size(); i++){
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), [&](int a, int b){
        return e->at(a) > e->at(b);
    });
    std::cout << "sorted by energy" << std::endl;
    std::vector<bool> used(e->size(), false);

    for (size_t j = 0; j < indices.size(); j++){
        int idx = indices[j];
        if (used[idx]) continue;

        for (size_t k = j + 1; k < indices.size(); k++){
            int idx2 = indices[k];
            if (used[idx2]) continue;
	    if (!t||!x||!y||!z||!e||!om) return;
            // Check time proximity
	   // if (idx >= t->size() || idx2 >= t->size() || idx >= x->size() || idx2 >= x->size() || idx >= y->size() ||
	//	idx2 >= y->size() || idx >= z->size() || idx2 >= z->size() || idx>= e->size() || idx2>=e->size()) continue;
            if (std::abs(t->at(idx) - t->at(idx2)) > max_dt) continue;
	    std::cout << "within time window" << std::endl;
            // Check sign change
	    if (idx >= x->size()|| idx >= y->size() || idx >= z->size() || idx2 >= x->size() || idx2 >= y->size() ||
		idx2 >= z->size()) return;
            bool oppx = (x->at(idx)*x->at(idx2) < 0);
            bool oppy = (y->at(idx)*y->at(idx2) < 0);
            bool oppz = (z->at(idx)*z->at(idx2) < 0);

            if (oppx || oppy || oppz){
                eout.push_back(e->at(idx));
                eout.push_back(e->at(idx2));
                tout.push_back(t->at(idx));
                tout.push_back(t->at(idx2));
                xout.push_back(x->at(idx));
                xout.push_back(x->at(idx2));
                yout.push_back(y->at(idx));
                yout.push_back(y->at(idx2));
                zout.push_back(z->at(idx));
                zout.push_back(z->at(idx2));
		omout.push_back(om->at(idx));
		omout.push_back(om->at(idx2));
                used[idx] = true;
                used[idx2] = true;
                break;
            }
	   
        }
    }
    std::cout << "New number of clusters: " << eout.size() << std::endl;
}


void angularDistribution(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
TH2D* h1, TH1D* h2, TH1D* h3){
    // Calculate distance between first and last cluster centres
    if (x.empty() || y.empty() || z.empty()) return;
    if (x.size() != y.size() || y.size() != z.size() || z.size() != x.size()) return;
    double dx = x.front() - x.back();
    double dy = y.front() - y.back();
    double dz = z.front() - z.back();
    // std::cout << "dx: " << dx << std::endl;
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    // Calculate angles
    double theta = acos(dz/r)*(180/TMath::Pi());
    double phi = atan2(dy, dx)*(180/TMath::Pi());
    // Fill histograms
    h1->Fill(theta, phi);
    h2->Fill(theta);
    h3->Fill(phi);
}


void checkDirection(const std::vector<double>& x, std::vector<double>& tdc, int p2n, int n2p, TH1D* h4, TH1D* h5){
    p2n = 0;
    n2p = 0;
    if (x.empty() || x.size() < 2) return;
    double dt = tdc.back() - tdc.front();
    // If x coord of cluster is negative, the cluster is on the Italian wall
    // If positive, then the cluster is on the French wall
    // If x1 < 0, x2 > 0: IT->FR
    // If x1 > 0, x1 < 0: FR->IT
    if (x.front() < 0 && x.back() > 0){
        n2p++;
        h4->Fill(dt);
    }
    else if (x.front() > 0 && x.back() < 0){
        p2n++;
        h5->Fill(dt);
    }
}

void ratioPlot(const std::vector<double>& x, const std::vector<double>& z, TH1D* h6, TH1D* h7){
    if (x.front() < 0 && x.back() > 0){
	h6->Fill(z.front());
	h6->Fill(z.back());
    }
    else if (x.front() > 0 && x.back() < 0){
	h7->Fill(z.front());
	h7->Fill(z.back());
    }
}

void ratioPlot2(const std::vector<double>& x, const std::vector<int>& om_number, TH1I* h8, TH1I* h9){
   if (x.front() < 0 && x.back() > 0){
	h8->Fill(om_number.front());
	h8->Fill(om_number.back());
   }
   else if (x.front() > 0 && x.back() < 0){
	h9->Fill(om_number.front());
        h9->Fill(om_number.back());
   }
} 

void caloAnalysis3(){
    gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
    TString parent = "/sps/nemo/scratch/Edinburgh_student/Katherine/CALO3/"; // path to folder containing individual run result folders
    TString prefix = "run-"; // prefix of every directory name
    TString filename = "clusteroutput.root"; // common file name
    TString treeName = "Result_tree";

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
    // Variables
    std::vector<double>* c_energy = nullptr;
    std::vector<std::vector<int>>* c_OMnums = nullptr;
    std::vector<double>* c_XC = nullptr;
    std::vector<double>* c_YC = nullptr;
    std::vector<double>* c_ZC = nullptr;
    std::vector<int>* c_centreOM = nullptr;
    std::vector<double>* c_centreTDC = nullptr;
    std::vector<int>* c_wallID = nullptr;
    chain->SetBranchAddress("c_energy", &c_energy);
    chain->SetBranchAddress("c_OMnums", &c_OMnums);
    chain->SetBranchAddress("c_XC", &c_XC);
    chain->SetBranchAddress("c_YC", &c_YC);
    chain->SetBranchAddress("c_ZC", &c_ZC);
    chain->SetBranchAddress("c_centreOM", &c_centreOM);
    chain->SetBranchAddress("c_centreTDC", &c_centreTDC);
    chain->SetBranchAddress("c_wallID", &c_wallID);
    // Histogram objects
    TH2D *h1 = new TH2D("h1", "2D Angular distribution", 16, 0, 180, 32, -180, 180);
    TH1D *h2 = new TH1D("h2", "#theta", 100, -180, 180);
    TH1D *h3 = new TH1D("h3", "#phi", 100, -180, 180);
    TH1D *h4 = new TH1D("h4", "IT->FR/FR->IT time difference", 100, 0, 50);
    TH1D *h5 = new TH1D("h5", "FR->IT time difference", 100, 0, 50);
    TH1D *h6 = new TH1D("h6", "IT->FR/FR->IT ratio", 259, -2000, 2000);
    TH1D *h7 = new TH1D("h7", "h7", 259, -2000, 2000);
    TH1I *h8 = new TH1I("h8", "IT->FR/FR->IT ratio", 712, 0, 712);
    TH1I *h9 = new TH1I("h9", "h9", 712, 0, 712);
    // Counter for IT->FR events vs FR->IT events
    int itfr = 0;
    int frit = 0;
    // New vectors
    std::vector<double> e_out;
    std::vector<double> t_out;
    std::vector<double> x_out;
    std::vector<double> y_out;
    std::vector<double> z_out;
    std::vector<int> om_out;
    Long64_t n = chain->GetEntries();
    // std::cout << "Number of entries: " << n << std::endl;
    for (Long64_t i = 0; i < n; i++){
        chain->GetEntry(i);
        std::cout << i << std::endl;
	std::cout << chain->GetCurrentFile()->GetName() << std::endl;
        // Apply cut on two highest energy clusters
        cutClusters(c_energy, c_centreTDC, c_XC, c_YC, c_ZC, c_centreOM, e_out, t_out, x_out, y_out, z_out, om_out, 20.0);
        // Sort by time
        sortallVectors(x_out, y_out, z_out, e_out, t_out, om_out);
        // for (size_t j = 0; j < y_out.size(); j++){
        //     std::cout << t_out.at(j) << std::endl;
        // }
        // Calculate angles
	std::cout << "cut and sorted clusters" << std::endl;
        angularDistribution(x_out, y_out, z_out, h1, h2, h3);
	std::cout << "calculated angles" << std::endl;
        checkDirection(x_out, t_out, 0, 0, h4, h5);
	ratioPlot(x_out, z_out, h6, h7);
	ratioPlot2(x_out, om_out, h8, h9);
    }
    std::cout << "Number of entries: " << n << std::endl;

    // Plotting
    TCanvas *c1 = new TCanvas("c1", "Angular distributions", 1600, 1200);
    TCanvas *c2 = new TCanvas("c2", "IT->FR/FR->IT time diff", 1600, 1200);
    TCanvas *c3 = new TCanvas("c3", "IT->FR/FR->IT ratio plot", 1600, 1200);
    TCanvas *c4 = new TCanvas("c4", "IT->FR/FR->IT ratio plot2", 1600, 1200);
    c1->cd();
    c1->Divide(2,2);
    c1->cd(1);
    h1->SetStats(0);
    h1->Draw("COLZ");
    gStyle->SetPalette(kLake);
    h1->GetXaxis()->SetTitle("Elevation angle #theta (degrees)");
    h1->GetYaxis()->SetTitle("Azimuthal angle #phi (degrees)");
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleOffset(0.95);
    h1->GetXaxis()->SetLabelSize(0.04);
    h1->GetYaxis()->SetLabelSize(0.04);
    c1->Modified();
    c1->Update();
    c1->cd(2);
    h2->SetStats(0);
    h2->SetMarkerStyle(2);
    h2->Draw();
    h2->GetXaxis()->SetTitle("Elevation angle #theta (degrees)");
    h2->GetYaxis()->SetTitle("Entries");
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->GetXaxis()->SetLabelSize(0.04);
    h2->GetYaxis()->SetLabelSize(0.04);
    c1->Modified();
    c1->Update();
    double meanTheta = h2->GetMean();
    int binmax = h2->GetMaximumBin();
    double maxTheta = h2->GetXaxis()->GetBinCenter(binmax);
    double maxThetaerr = h2->GetBinError(binmax);
    std::cout << "Zenith angle mode: " << maxTheta << " +/- " << maxThetaerr << " degrees" << std::endl;
    std::cout << "Mean zenith angle: " << meanTheta << " degrees" << std::endl;
    TF1 *fit = new TF1("fit", "gaus(0) + gaus(3)", 0, 180);
    fit->SetParameter(0, 20); // amp 1
    fit->SetParameter(1, 45); // mean 1
    fit->SetParameter(2, 40); // sigma 1
    fit->SetParameter(3, 20); // amp 2
    fit->SetParameter(4, 135); // mean 2
    fit->SetParameter(5, 4); // sigma 2
    h2->Fit(fit, "R");
    fit->SetLineColor(kBlue);
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
    c1->Modified();
    c1->Update();
    double meanPhi = h3->GetMean();
    int binmax2 = h3->GetMaximumBin();
    double maxPhi = h3->GetXaxis()->GetBinCenter(binmax2);
    double maxPhierr = h3->GetBinError(binmax2);
    std::cout << "Azimuth angle mode: " << maxPhi << "+/-" << maxPhierr << " degrees" << std::endl;
    std::cout << "Mean azimuth angle: " << meanPhi << " degrees" << std::endl;
    c1->SaveAs("caloangdistributionfitted.png");
    c2->cd();
    h4->SetStats(0);
    h4->Draw();
    h4->GetXaxis()->SetTitle("Time difference #Delta t (ns)");
    h4->GetYaxis()->SetTitle("Entries");
    h5->SetStats(0);
    h5->Draw("SAME");
    h5->SetLineColor(kRed);
    auto legend2 = new TLegend(0.7, 0.75, 0.9, 0.85);
    legend2->SetTextSize(0.02);
    legend2->AddEntry(h4, "IT->FR events", "l");
    legend2->AddEntry(h5, "FR->IT events", "l");
    legend2->Draw();
    c2->SaveAs("ratiotimediff.png");
    c3->cd();
    TRatioPlot *rp = new TRatioPlot(h6, h7, "divsym");
    h6->SetStats(0);
    double int1 = h6->Integral();
    // h6->Scale(1/h6->Integral());
   // h6->SetMarkerStyle(21);
    h6->SetMarkerColor(kRed);
    h6->SetLineColor(kRed);
    h6->GetXaxis()->SetTitle("z position");
   // h7->SetMarkerStyle(21);
    h7->SetMarkerColor(kBlack);
    // h7->Scale(1/h7->Integral());
    double int2 = h7->Integral();
    h7->SetLineColor(kBlack);
    rp->SetH1DrawOpt("PE");
    rp->SetH2DrawOpt("PE");
    rp->Draw();
    rp->GetLowerRefYaxis()->SetTitle("Ratio");
    rp->GetUpperRefYaxis()->SetTitle("Entries");
    rp->GetUpperPad()->cd();
    rp->GetLowYaxis()->SetNdivisions(505);
    TLegend *leg = new TLegend(0.7, 0.75, 1.0, 0.85);
    leg->AddEntry("h6", "IT->FR events", "l");
    leg->AddEntry("h7", "FR->IT events", "le");
    leg->Draw();
    c3->Update();
    c3->SaveAs("clusterdirectionalityrationorm.png");
    double ratio = int1/int2; // raw ratio of It->FR/FR->IT
    std::cout << "Ratio = " << ratio << std::endl;
    c4->cd();
    TRatioPlot *rp2 = new TRatioPlot(h8, h9, "divsym");
    h8->SetStats(0);
    h8->SetMarkerColor(kRed); 
    h8->SetLineColor(kRed);
    h8->GetXaxis()->SetTitle("OM number");
    h9->SetLineColor(kBlack);
    h9->SetMarkerColor(kBlack);
    rp2->SetH1DrawOpt("PE");
    rp2->SetH2DrawOpt("PE");
    rp2->Draw();
    rp2->GetLowerRefYaxis()->SetTitle("Ratio");
    rp2->GetUpperRefYaxis()->SetTitle("Entries");
    rp2->GetUpperPad()->cd();
    rp2->GetLowYaxis()->SetNdivisions(505);
    TLegend *lego = new TLegend(0.7, 0.75, 1.0, 0.85);
    lego->AddEntry("h8", "IT->FR events", "l");
    lego->AddEntry("h9", "FR->IT events", "le");
    lego->Draw();
    c4->Update();
    c4->SaveAs("OMratioplot.png");
    // Output file
    TFile *f = new TFile("calohistograms.root", "RECREATE");
    h1->Write();
    h2->Write();
    h3->Write();
    f->Close();
}
