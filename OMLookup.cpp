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
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TList.h>

struct coordinate{
    double x;
    double y;
    double z;
    bool nonnullCoord = false; // false by default, tells the lookup script whether value non-null
};



void OMLookup(){
    TChain *chain = new TChain("Result_tree");
    chain->Add("/sps/nemo/scratch/Edinburgh_student/Katherine/data/reconstructed_runs_old/run-*.root");
    // Variables (om_number, gamma_om_x,y,z)
    std::vector<int>* omnumber = nullptr;
    std::vector<double>* gamma_om_x = nullptr;
    std::vector<double>* gamma_om_y = nullptr;
    std::vector<double>* gamma_om_z = nullptr;
    // New variables for new tree
    std::vector<int> om_number;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    // Branches
    chain->SetBranchAddress("om_number", &omnumber);
    chain->SetBranchAddress("gamma_om_x", &gamma_om_x);
    chain->SetBranchAddress("gamma_om_y", &gamma_om_y);
    chain->SetBranchAddress("gamma_om_z", &gamma_om_z);
    Long64_t n = chain->GetEntries();
    std::map<int, coordinate> coordinateMap;
    // Event loop
    for (Long64_t i = 0; i < n; i++){
        chain->GetEntry(i);
        om_number.clear();
        x.clear();
        y.clear();
        z.clear();
        if (!omnumber || !gamma_om_x || !gamma_om_y || !gamma_om_z) continue;
        for (int j = 0; j < gamma_om_x->size(); j++){
            auto [it, in] = coordinateMap.emplace(0, coordinate{-450.5, -2460.5, -1554});
            int current_om = omnumber->at(j);
                if (coordinateMap.find(current_om) == coordinateMap.end()){
                coordinate c;
                if (gamma_om_x->at(j) != -9999.99 || gamma_om_y->at(j) != -9999.99 || gamma_om_z->at(j) != -9999.99){
                    c.nonnullCoord = true;
                    c.x = gamma_om_x->at(j);
                    c.y = gamma_om_y->at(j);
                    c.z = gamma_om_z->at(j);
                    coordinateMap[current_om] = c;
                }
                else{
                    c.nonnullCoord = false;
                }
                
            }
            
        }
        
    }
    // Write to a .csv file
    std::ofstream fout("/sps/nemo/scratch/Edinburgh_student/Katherine/OMGeometry.csv");
    fout << "# om_number \t x \t y \t z\n";
    for (const auto& [om, c] : coordinateMap){
        double X = c.x;
        double Y = c.y;
        double Z = c.z;
        fout << om << "\t" << X << "\t" << Y << "\t" << Z << "\n";
    }
    fout.close();
}
