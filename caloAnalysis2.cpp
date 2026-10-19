#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <queue>
#include <map>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TObject.h>
#include <TROOT.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TPolyLine3D.h>
#include <Fit/Fitter.h>
#include <Math/Vector3D.h>

// void sortTDC(std::vector<double>* vec, std::vector<double>* tdc){
    
//     // std::sort(vec->begin(), vec->end())
// }

// void ITFRratio(std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, std::vector<double>* tdc){
//     // Order events by increasing tdc
//     std::vector<int> idx;
//     std::iota(idx.begin(), idx.end(), 0);
//     std::sort(idx.begin(), idx.end(), [&](int a, int b){
//         return tdc[a] < tdc[b];
//     });
//     for (size_t j = 0; j < x->size(); j++){
//         std::cout << tdc->at(j) << std::endl;
//     }
    
// }
void printMatrices(const std::vector<std::vector<double>>& matrix, int eventNumber, const char* branch){
    const size_t l = matrix.size();
    // Column width
    constexpr int width = 12;
    // Column headers
    std::cout << std::setw(width) << "Event " + std::to_string(eventNumber) + " "+ branch + ": |"; //spacer
    for (int o = 0; o < l; o++){
        std::cout << std::setw(width) << ("Cluster " + std::to_string(o));
    }
    std::cout << "\n";
    // Separator
    std::cout << std::string(width, '-') << "-+-" << std::string(l*width, '-') << "\n";
    // Print row labels, data and separators
    for (int p = 0; p < l; p++){
        // Header
        std::cout << std::setw(width) << ("Cluster " + std::to_string(p)) << " |";

        // Matrix elements
        for (int q = 0; q < l; q++){
            std::cout << std::setw(width) << std::fixed << std::setprecision(2) << matrix[p][q];
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

void makeweirdMatrices(const std::vector<std::vector<int>>* source, std::vector<std::vector<double>>& matrix){
    int m = source->size();
    matrix.clear();
    // matrix.assign(m, std::vector<double>(m, 0.0));
    matrix.resize(m, std::vector<double>(m, 0.0));
    for (int i = 0; i < m; i++){
        for (int j = 0; j < m; j++){
            
            size_t sizeI = static_cast<int>(source[i].size());
            size_t sizeJ = static_cast<int>(source[j].size());

            std::cout << sizeI << std::endl;
            std::cout << source[i].size() << std::endl;
            matrix[i][j] = static_cast<double>(sizeJ - sizeI);
            // matrix[i][j] = static_cast<double>(std::abs(static_cast<int>(source[i].size()) - static_cast<int>(source[j].size())));
            }
            
        }
    }


// void makeweirdMatrices(const std::vector<std::vector<int>>* source, std::vector<std::vector<double>>& matrix){
//     int rows = source->size();
//     matrix.clear();
//     if (rows < 2) return;
//     matrix.resize(rows - 1, std::vector<double>(1, 0.0));
//     for (int i = 0; i < rows - 1; i++){
//         int sizeI = static_cast<int>(source[i].size());
//         int sizeJ = static_cast<int>(source[i+1].size());
//         matrix[i][0] = static_cast<double>(sizeJ - sizeI);
//     }
// }

void clusterCSV(int evnum, std::vector<double>* x, std::vector<double>* y, std::vector<double>* z,
std::vector<double>* energy, std::vector<int> clusterSize, std::vector<double>* tdc,TString folder){
    TString name = Form("%s/clusterdata_%d.csv", folder.Data(),evnum);
    std::ofstream fout(name);
    fout << "event_number \t XC \t YC \t ZC \t Cluster energy \t Cluster size \t Centre tdc \n";
    for (size_t k = 0; k < x->size(); k++){
        fout <<  x->at(k) << "," << y->at(k) << "," << z->at(k) << "," << energy->at(k) << "," << clusterSize.at(k) << "," << tdc->at(k) <<  "," << evnum << "\n";
    }
    fout.close();
}

void makeMatrices(std::vector<double>* vec, std::vector<std::vector<double>>& matrix){
    // Matrix vector
    int m = vec->size();
    matrix.assign(m, std::vector<double>(m));
    for (size_t j = 0; j < vec->size(); j++){
        for (size_t k = j; k < vec->size(); k++){
            double diff = std::abs(vec->at(j) - vec->at(k));
            matrix[j][k] = diff;
            matrix[k][j] = diff;

        }
    }
}

void distance(const std::vector<double>* x, const std::vector<double>* y, const std::vector<double>* z, std::vector<double>& r){
    for (size_t j = 0; j < x->size(); j++){
        double X = x->at(j);
        double Y = y->at(j);
        double Z = z->at(j);
        double rc = sqrt(X*X + Y*Y + Z*Z);
        r.push_back(rc);
    }
}

void makedistanceMatrix(std::vector<double>& vec, std::vector<std::vector<double>>& matrix){
    // Matrix vector
    int m = vec.size();
    matrix.assign(m, std::vector<double>(m, 0.0));
    for (size_t j = 0; j < vec.size(); j++){
        for (size_t k = j; k < vec.size(); k++){
            double diff = std::abs(vec.at(j) - vec.at(k));
            matrix[j][k] = diff;
            matrix[k][j] = diff;

        }
    }
}



void caloAnalysis2(){
    gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
    TString parent = "/sps/nemo/scratch/Edinburgh_student/Katherine/CALO2/"; // path to folder containing individual run result folders
    TString prefix = "run-"; // prefix of every directory name
    TString filename = "clusteroutput.root"; // common file name
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

    // Variables
    std::vector<double>* c_energy = nullptr;
    std::vector<double>* c_centreTDC = nullptr;
    std::vector<std::vector<int>>* c_OMnums = nullptr;
    int event_number;
    std::vector<double>* c_XC = nullptr;
    std::vector<double>* c_YC = nullptr;
    std::vector<double>* c_ZC = nullptr;
    // Branch addresses
    chain->SetBranchAddress("c_energy", &c_energy);
    chain->SetBranchAddress("c_centreTDC", &c_centreTDC);
    chain->SetBranchAddress("c_OMnums", &c_OMnums);
    chain->SetBranchAddress("event_number", &event_number);
    chain->SetBranchAddress("c_XC", &c_XC);
    chain->SetBranchAddress("c_YC", &c_YC);
    chain->SetBranchAddress("c_ZC", &c_ZC);
    // Instantiate matrices
    std::vector<std::vector<double>> energyMatrix;
    std::vector<std::vector<double>> tdcMatrix;
    std::vector<std::vector<double>> sizeMatrix;
    std::vector<std::vector<double>> distMatrix;
    Long64_t n = chain->GetEntries();
    TString branch = "energy";
    TString branch2 = "TDC";
    TString branch3 = "cluster size";
    TString branch4 = "distance";
    TString skib = "skibidi";
    for (Long64_t i = 0; i < n; i++){
        //Long64_t entryNumber = chain->GetEntryNumber(i);
        //std::cout << skib << std::endl;
        TString dirPath = gSystem->DirName(chain->GetCurrentFile()->GetName());
        std::cout << "FILE: " << dirPath << std::endl;
        // loop over events
        chain->GetEntry(i);
        makeMatrices(c_energy, energyMatrix);
        printMatrices(energyMatrix, event_number, branch);
        makeMatrices(c_centreTDC, tdcMatrix);
        printMatrices(tdcMatrix, event_number, branch2);
        // makeweirdMatrices(c_OMnums, sizeMatrix);
        // printMatrices(sizeMatrix, event_number, branch3);
        std::vector<double> r;
        distance(c_XC, c_YC, c_ZC, r);
        makedistanceMatrix(r, distMatrix);
        printMatrices(distMatrix, event_number, branch4);
        size_t nClust = c_energy->size();
        std::cout << "Found " << nClust << "Clusters" << std::endl;
        // we want to write to the CSV the number of hits / OM happenings in each event,
        // so c_OMnums needs to be replaced with an array that looks like
        // [N_hits_cluster1,N_hits_cluster2, # hits in cluster 3]
        std::vector<int> ClusterSize(nClust);
        for (size_t clusterCount = 0; clusterCount < nClust; clusterCount++){
            std::cout << clusterCount;
            size_t hitNum = c_OMnums->at(clusterCount).size();
            std::cout << " : " << hitNum << std::endl;
            ClusterSize[clusterCount] = hitNum;
        }
        clusterCSV(event_number, c_XC, c_YC, c_ZC, c_energy, ClusterSize, c_centreTDC,dirPath);
        //sortTDC(c_XC, c_centreTDC);
        // ITFRratio(c_XC, c_YC, c_ZC, c_centreTDC);
        //tree->Print();
        //tree->Scan("c_OMnums");
        // for (const auto& row : *c_OMnums){
            // std::cout << row.size() << std::endl;
        // }
    }
    // energyTable(tree);
}
