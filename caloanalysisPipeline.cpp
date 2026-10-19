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

const int N_rows = 13; // number of rows of OMs in MWall

struct OMHit {
    int omnum; // ID of the OM hit
    double energy; // energy of the hit
    double tdc; // calo_tdc of hit
    double x;
    double y;
    double z;
    bool isHit = false;

    // Functions to get coordinates from the OM number
    int U() const{return omnum % N_rows;} // "X coordinate" of OM
    int V() const{return omnum / N_rows;} // "Y coordinate" of OM

};

struct OMCluster{
    double totEnergy = 0; // total energy of cluster
    double U_COM = 0; // COM of cluster U coord
    double V_COM = 0; // COM of cluster V coord
    std::vector<int> OMIDs;
    std::vector<double> OM_tdc;
    int centreOM; // ID of the OM closest to the centre of the cluster
    double centreOMTDC; // TDC value corresponding to centre
    int wallID; // designates which wall a hit occurs on (0 for Italy, 1 for France)
    double XC = 0.0;
    double YC = 0.0;
    double ZC = 0.0;
};

struct coordinate{
    double x;
    double y;
    double z;
};

// bool neighbourID(const OMHit& h1, const OMHit& h2){
//     // Checks if two OMs are neighbours
//     int dU = std::abs(h1.U() - h2.U());
//     int dV = std::abs(h1.V() - h2.V());
//     return (dU <= 1 && dV <= 1); // true if their edges or corners touch
// }

double distSquared(const OMHit& h1, const OMHit& h2){
    double dx = h1.x - h2.x;
    double dy = h1.y - h2.y;
    double dz = h1.z - h2.z;
    return (dx*dx + dy*dy + dz*dz); 
}


std::vector<OMCluster> FindClusters(std::vector<OMHit>& hits, double eThresh, double maxTime, double neighbourThresh){
    std::vector<OMCluster> clusters;
    // Mapping OM hit IDs to an index in the hits vector
    std::map<int, size_t> IDidxmap;
    for (size_t l = 0; l < hits.size(); l++){
        IDidxmap[hits[l].omnum] = l;
    }

    for (size_t m = 0; m < hits.size(); m++){
        if (hits[m].energy < eThresh || hits[m].isHit) continue;

        OMCluster cluster;
        std::queue<size_t> q;
        q.push(m);
        hits[m].isHit = true;

        double EsumU = 0;
        double EsumV = 0;

        while (!q.empty()){
            size_t current = q.front();
            q.pop();

            OMHit& currentHit = hits[current];
            // if (abs(currentHit.tdc) >= maxTime) continue;
            cluster.OMIDs.push_back(currentHit.omnum);
            cluster.OM_tdc.push_back(currentHit.tdc);
            cluster.totEnergy += currentHit.energy;

            // Weighting cluster COM
            EsumU += currentHit.energy * currentHit.U();
            EsumV += currentHit.energy * currentHit.V();

            // Finding neighbours
            int currentU = currentHit.U();
            int currentV = currentHit.V();
            for (int dU = -1; dU <= 1; dU++){
                for (int dV = -1; dV <= 1; dV++){
                    if (dU == 0 && dV == 0) continue; // skip current OM
                    int neighbourU = currentU + dU;
                    int neighbourV = currentV + dV;
                    // Calculate the OM numbers of neighbouring OMs
                    int neighbourOM = neighbourV*N_rows + neighbourU;
                    // Check if neighbouring OMs were hit
                    auto it = IDidxmap.find(neighbourOM);
                    if (it != IDidxmap.end()){
                        size_t neighbouridx = it->second;
                        if (!hits[neighbouridx].isHit && hits[neighbouridx].energy >= neighbourThresh){
                            double dt = std::abs(currentHit.tdc - hits[neighbouridx].tdc);
                            if (dt <= maxTime){
                                hits[neighbouridx].isHit = true;
                                q.push(neighbouridx);
                            }
                            // hits[neighbouridx].isHit = true;
                            // q.push(neighbouridx);
                        }
                    }
                }
            }
        }
        if (cluster.totEnergy > 0){
            cluster.U_COM = EsumU/cluster.totEnergy;
            cluster.V_COM = EsumV/cluster.totEnergy;
            // Find the OM number closest to the centre of the cluster
            int centreU = std::round(cluster.U_COM);
            int centreV = std::round(cluster.V_COM);
            cluster.centreOM = (centreV*N_rows) + centreU;
            if (cluster.centreOM < 260){
                cluster.wallID = 0;
            }
            else if (cluster.centreOM > 260){
                cluster.wallID = 1;
            }
            // Look up time of centre hit
            auto cit = IDidxmap.find(cluster.centreOM);
            if (cit != IDidxmap.end()) {
                cluster.centreOMTDC = hits[cit->second].tdc;
            }
            // cluster.centreOMTDC = cluster.OM_tdc[cluster.centreOM];
        }
        clusters.push_back(cluster);
    }
    return clusters;  
}

struct omcoordinate{
    double OMx;
    double OMy;
    double OMz;
};


void clusterAnalysis(TTree* cal_tree){
    gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
    // Open .csv file containing coordinates of the centres of each OM (see omgeometry.cpp)
    std::ifstream inFile("/sps/nemo/scratch/Edinburgh_student/Katherine/OMGeometry.csv");
    std::map<int, coordinate> OMLookup;
    std::string line;
    std::getline(inFile, line);
    int om;
    double x, y, z;
    while (inFile >> om >> x >> y >> z){
        coordinate c = {x, y, z};
        OMLookup[om] = c;
    }
    // Checking map has loaded correctly
    // for (const auto& [id, p] : OMLookup){
    //     std::cout << id << ", coords: " << p.x << std::endl;
    // }
    inFile.close();
    // Variables
    int OM_count;
    std::vector<double>* energy = nullptr;
    std::vector<int>* om_number = nullptr;
    int event_number;
    std::vector<double>* cal_calo_tdc = nullptr;
    std::vector<double>* gamma_om_x = nullptr; // x, y, z coordinates for OMs where a gamma has interacted
    std::vector<double>* gamma_om_y = nullptr; // doesn't have to be in gamma veto 
    std::vector<double>* gamma_om_z = nullptr;
    // sortVectors(calo_tdc, om_number, calo_tdc);
    // Branches
    cal_tree->SetBranchAddress("OM_count", &OM_count);
    cal_tree->SetBranchAddress("energy", &energy);
    cal_tree->SetBranchAddress("om_number", &om_number);
    cal_tree->SetBranchAddress("event_number", &event_number);
    cal_tree->SetBranchAddress("cal_calo_tdc", &cal_calo_tdc);
    cal_tree->SetBranchAddress("gamma_om_x", &gamma_om_x);
    cal_tree->SetBranchAddress("gamma_om_y", &gamma_om_y);
    cal_tree->SetBranchAddress("gamma_om_z", &gamma_om_z);
    // int len = 0;
    // Cluster output files and output trees
    // File for just events with multiple clusters
    TFile *fOut = new TFile("clusteroutput.root", "RECREATE");
    TTree *outTree = cal_tree->CloneTree(0);
    // Cluster variables
    std::vector<double> c_energy; // total energy of cluster
    std::vector<double> c_U; // x coords of cluster
    std::vector<double> c_V; // y coords of cluster
    std::vector<std::vector<int>> c_OMnums; // OM numbers of hits in cluster
    std::vector<std::vector<double>> c_OMtdc; // tdc values of all hits in cluster
    std::vector<int> c_centreOM; // OM numbers of OMs closest to cluster centre
    std::vector<double> c_centreTDC; // TDC values for OMs at cluster centres
    std::vector<int> c_wallID; // wall IDs for clusters
    std::vector<double> c_XC; // x coordinates of the centre of the cluster
    std::vector<double> c_YC; // y coordinates
    std::vector<double> c_ZC; // z coordinates
    double theta = 0;
    double phi = 0;
    // Setting output branches
    TBranch *b1 = outTree->Branch("c_energy", &c_energy);
    TBranch *b2 = outTree->Branch("c_U", &c_U);
    TBranch *b3 = outTree->Branch("c_V", &c_V);
    TBranch *b4 = outTree->Branch("c_OMnums", &c_OMnums);
    TBranch *b5 = outTree->Branch("c_centreOM", &c_centreOM);
    TBranch *b6 = outTree->Branch("c_centreTDC", &c_centreTDC);
    TBranch *b7 = outTree->Branch("c_wallID", &c_wallID);
    TBranch *b8 = outTree->Branch("c_XC", &c_XC);
    TBranch *b9 = outTree->Branch("c_YC", &c_YC);
    TBranch *b10 = outTree->Branch("c_ZC", &c_ZC);
    TBranch *b11 = outTree->Branch("event_number", &event_number);
    // TBranch *b11 = outTree->Branch("theta", &theta);
    // TBranch *b12 = outTree->Branch("phi", &phi);
    // Cluster vectors 
    std::vector<std::vector<OMCluster>> multiClusters; // Vector to store multi-cluster events
    std::vector<std::vector<OMCluster>> allClusters; // Vector to store all clusters
    std::vector<std::vector<OMCluster>> ITFRClusters; // Events going from IT->FR
    std::vector<std::vector<OMCluster>> FRITClusters; // Events going from FR->IT
    // Histogram objects
    TH1D *H1 = new TH1D("H1", "Event total energy", 100, 0, 100);
    TH1D *H2 = new TH1D("H2", "Cluster total energy", 100, 0, 100);
    TH2I *H3 = new TH2I("H3", "Energies of individual OMs", 89, 0, 712, 100, 0, 100);
    TH2D *H4 = new TH2D("H4", "Cluster timing vs energy", 50, 0, 6000, 50, 0, 100);
    TH2D *H5 = new TH2D("H5", "Cluster timing vs size", 50, 0, 6000, 20, 0, 20);
    TH1D *H6 = new TH1D("H6", "Cluster timing distribution", 100, 0, 6000);
    // Event loop
    Long64_t n = cal_tree->GetEntries();
    for (Long64_t i = 0; i < n; i++){   
        cal_tree->GetEntry(i);
        
        // Clear previous entry
        c_energy.clear();
        c_U.clear();
        c_V.clear();
        c_OMnums.clear();
        c_OMtdc.clear();
        c_centreOM.clear();
        c_centreTDC.clear();
        c_wallID.clear();
        c_XC.clear();
        c_YC.clear();
        c_ZC.clear();
        
        // Convert data to OMHit struct format
        std::vector<OMHit> Hits;
        size_t noHits = om_number->size(); // number of hits in event
        Hits.reserve(noHits);
        for (size_t g = 0; g < noHits; g++){
            auto it = OMLookup.find(om_number->at(g));
            if (it != OMLookup.end()){
                Hits.push_back({om_number->at(g), energy->at(g), cal_calo_tdc->at(g), it->second.x, it->second.y, it->second.z, false});
            }
            
            // if (om_number->at(g) > 520) continue; // use MWall only (needed if not using time calibration)
        }
        // try{
        //     auto [p1, p2] = furthestHits(Hits);
        //     std::cout << i << ": Furthest hits: (" << p1.U() << ", " << p1.V() << ") and (" << p2.U() << ", " << p2.V() << ")" << std::endl;
        // } catch (const std::exception& e){
        //     std::cerr << e.what() << std::endl;
        // }
        // Mock hits for the two corners of the calorimeter
        OMHit h1 = {0, 0.0, 0.0, -450.5, -2460.5, -1554, true};
        OMHit h2 = {519, 0.0, 0.0, 450.5, 2460.5, 1554, true};
        // Find distance between the two corners
        double maxDistSq = distSquared(h1, h2);
        double maxDist = sqrt(maxDistSq);
        // std::cout << maxDist << std::endl;
        // Find time it takes to cross calorimeter
        double c = 300.0; // speed of light in mm/ns
        double mTime = maxDist/(0.99*c);
        // std::cout << mTime << " ns" << std::endl;
        // Find clusters using clustering function, for now set threshold energy as 0.0001 (just want non-zero)
        
        auto eventClusters = FindClusters(Hits, 1, 25, 0.00001);
        size_t nClust = eventClusters.size();
        std::cout << "L291: nClust THIS IS NUMBER OF CLUSTERS " << nClust << std::endl;
        // Fill output variables
        for (auto& c : eventClusters){
            // Add the energies and centre of mass
            c_energy.push_back(c.totEnergy);
            c_U.push_back(c.U_COM);
            c_V.push_back(c.V_COM);
            // add to c_OMnums, the array that contains 
            std::cout << "NUMBER OF OMs (and therefore hits?)  " << c.OMIDs.size() << std::endl;
            
            c_OMnums.push_back(c.OMIDs);
            // for (const auto& row : c_OMnums){
            //     for (int r : row){
            //         std::cout << "OM added: " << r << std::endl;
            //     }
            // }
            c_OMtdc.push_back(c.OM_tdc);
            // for (const auto& drow : c_OMtdc){
            //     for (double p : drow){
            //         std::cout << "TDC added: " << p << std::endl;
            //     }
            // }
            c_centreOM.push_back(c.centreOM);
            c_centreTDC.push_back(c.centreOMTDC);
            H5->Fill(c.centreOMTDC, c.OMIDs.size());
            c_wallID.push_back(c.wallID);
            // std::cout << i << ": Centre: " << c.centreOM << std::endl;
            // std::cout << i << ": Energy: " << c.totEnergy << std::endl;
            // std::cout << i << ": Number of OMs in cluster: " << c.OMIDs.size() << std::endl;
            // std::cout << i << ": Centre tdc: " << c.centreOMTDC << std::endl;
            

            auto a = OMLookup.find(c.centreOM);
            if (a != OMLookup.end()){
                
                c.XC = a->second.x;
                c.YC = a->second.y;
                c.ZC = a->second.z;
                c_XC.push_back(c.XC);
                c_YC.push_back(c.YC);
                c_ZC.push_back(c.ZC);
                // std::cout << i << ": " << c.centreOM << ": Energy: " << c.totEnergy << std::endl;
                // std::cout << i << ": " << c.centreOM << ": x: " << c.XC << std::endl;
                // std::cout << c.centreOM << ": y: " << c.YC << std::endl;
                // std::cout << c.centreOM << ": z: " << c.ZC << std::endl;
                
            }
            // Filling histograms
            double totE = 0.0;
            for (double w : *energy){
                totE += w;
                H1->Fill(totE);
            }
            for (size_t v = 0; v < c_energy.size(); v++){
                H2->Fill(c_energy.at(v));
            }
            for (size_t u = 0; u < om_number->size(); u++){
                H3->Fill(om_number->at(u), energy->at(u));
            }
            // size_t min_size = std::min(c_OMnums.size(), c_OMtdc.size());
            // for (size_t f = 0; f < c_OMnums.size(); f++){
            //     const auto& row1 = c_OMnums[f];
            //     const auto& row2 = c_OMtdc[f];
            //     for (int val : row1){std::cout << i << ": OM: " <<  val << std::endl;}
            //     for (double val2 : row2){std::cout << i <<  ": TDC: " <<  val2 << std::endl;}
            // }
            // for (const auto& row : c_OMtdc){
            //     for (double d : row){
            //         H4->Fill(d);
            //     }
            // }
            for (size_t d = 0; d < c_energy.size(); d++){
                H4->Fill(c_centreTDC.at(d), c_energy.at(d));
		H6->Fill(c_centreTDC.at(d));
            }
            // for (size_t ir = 0; ir < c_OMnums.size(); ir++){
            //     double cS = static_cast<double>(c_OMnums[ir].size());
            //     // std::cout << "Cluster size: " << cS << std::endl;
            //     H5->Fill(c_centreTDC.at(ir), cS);
            // }
            //for (size_t er = 0; er < c_OMnums.size(); er++){
               // std::cout << "NUMBER OF CLUSTERS?" << c_OMnums.size() << std::endl;
            //}
        
        }
        std::sort(eventClusters.begin(), eventClusters.end(), [](const OMCluster& a, const OMCluster& b){
            return a.centreOMTDC < b.centreOMTDC;
        });
        double noCluster = eventClusters.size(); // Number of clusters in an event
        std::cout << "Event " << i << " contains " << noCluster << " clusters." << std::endl;
        // Find direction
        // double dx = eventClusters.front().XC - eventClusters.back().XC;
        // double dy = eventClusters.front().YC - eventClusters.back().YC;
        // double dz = eventClusters.front().ZC - eventClusters.back().ZC;
        // double r = sqrt(dx*dx + dy*dy + dz*dz);
        // theta = acos(dz/r)*(180/TMath::Pi());
        // phi = atan2(dy, dx)*(180/TMath::Pi());
        // std::cout << theta << std::endl;
        
        // Fill tree - only keep events where there's more than 1 cluster
        fOut->cd();
        if (noCluster > 1){
            outTree->Fill();
        }
        
        // Find events where IT->FR and FR->IT
        if (eventClusters.front().wallID == 0 && eventClusters.back().wallID == 1){
            ITFRClusters.push_back(eventClusters);
        }
        else if (eventClusters.front().wallID == 1 && eventClusters.back().wallID == 0){
            FRITClusters.push_back(eventClusters);
        }
        allClusters.push_back(eventClusters);
        if (noCluster > 1){
            multiClusters.push_back(eventClusters);
        }
    }
    // Find fraction of multi-cluster events compared to all events
    double size1 = multiClusters.size();
    double size2 = allClusters.size();
    double frac = size1/size2;
    // Find ratio of IT->FR events to FR->IT events
    double sizeIT = ITFRClusters.size();
    double sizeFR = FRITClusters.size();
    double ratioITFR = sizeIT/sizeFR; // not sure how to represent this graphically
    std::cout << "Fraction of multi-cluster events: " << frac << std::endl;
    std::cout << "Ratio of IT->FR events to FR->IT events: " << ratioITFR << std::endl;
    TCanvas *c1 = new TCanvas("c1", "Event energy vs cluster energy", 1600, 1200);
    c1->Divide(2, 3);
    c1->cd(1);
    H1->SetStats(0);
    H1->Draw();
    H1->GetXaxis()->SetTitle("Event total energy (MeV)");
    H1->GetYaxis()->SetTitle("Entries");
    c1->cd(2);
    H2->SetStats(0);
    H2->Draw();
    H2->GetXaxis()->SetTitle("Cluster total energy (MeV)");
    H2->GetYaxis()->SetTitle("Entries");
    c1->cd(3);
    gStyle->SetPalette(kCool);
    H3->SetStats(0);
    H3->Draw("COLZ");
    H3->GetXaxis()->SetTitle("om_number");
    H3->GetYaxis()->SetTitle("Energy of hit (MeV)");
    c1->cd(4);
    H4->SetStats(0);    
    H4->Draw("COLZ");
    H4->GetXaxis()->SetTitle("Cluster timing (ns)");
    H4->GetYaxis()->SetTitle("Cluster energy (MeV)");
    H5->GetXaxis()->SetRangeUser(0, 4000);
    c1->cd(5);
    H5->SetStats(0);
    H5->Draw("COLZ");
    H5->GetXaxis()->SetTitle("Cluster timing (ns)");
    H5->GetYaxis()->SetTitle("Cluster size (number of OMs)");
    H5->GetXaxis()->SetRangeUser(0, 4000);
    c1->cd(6);
    H6->SetStats(0);
    H6->Draw();
    H6->GetXaxis()->SetTitle("Cluster timing (ns)");
    H6->GetYaxis()->SetTitle("Entries");
    H6->GetXaxis()->SetRangeUser(0, 4000);
    c1->SaveAs("clusterprelim.png");
    delete c1;
    TTree* finalTree = outTree->CopyTree("c_centreTDC <= 500");
    fOut->cd();
    if (finalTree->GetEntries() == 0){
        finalTree->Write(nullptr, TObject::kOverwrite);
    }
    else{
        finalTree->Write();
    }
    
    // fOut->Write();
    fOut->Close();
    // f->Close();
    

}


// void omgeometry(TTree* cal_tree){
//     // Variables (om_number, gamma_om_x,y,z)
//     std::vector<int>* omnumber = 0;
//     std::vector<double>* gamma_om_x = nullptr;
//     std::vector<double>* gamma_om_y = nullptr;
//     std::vector<double>* gamma_om_z = nullptr;
//     // New variables for new tree
//     std::vector<int> om_number;
//     std::vector<double> x;
//     std::vector<double> y;
//     std::vector<double> z;
//     // Branches
//     cal_tree->SetBranchAddress("om_number", &omnumber);
//     cal_tree->SetBranchAddress("gamma_om_x", &gamma_om_x);
//     cal_tree->SetBranchAddress("gamma_om_y", &gamma_om_y);
//     cal_tree->SetBranchAddress("gamma_om_z", &gamma_om_z);
//     Long64_t n = cal_tree->GetEntries();
//     std::map<int, coordinate> coordinateMap;
//     // Event loop
//     for (Long64_t i = 0; i < n; i++){
//         cal_tree->GetEntry(i);
//         om_number.clear();
//         x.clear();
//         y.clear();
//         z.clear();
//         if (!omnumber || !gamma_om_x || !gamma_om_y || !gamma_om_z) continue;
//         for (int j = 0; j < gamma_om_x->size(); j++){
//             int current_om = omnumber->at(j);
//             if (coordinateMap.find(current_om) == coordinateMap.end()){
//                 coordinate c;
//                 c.x = gamma_om_x->at(j);
//                 c.y = gamma_om_y->at(j);
//                 c.z = gamma_om_z->at(j);
//                 coordinateMap[current_om] = c;
//             }
//         }
//     }
//     // file->Close();
//     // Write to text file
//     std::ofstream fout("OMgeometry.txt");
//     fout << "# om_number \t x \t y \t z\n";
//     for (const auto& [om, c] : coordinateMap){
//         double X = c.x;
//         double Y = c.y;
//         double Z = c.z;
//         fout << om << "\t" << X << "\t" << Y << "\t" << Z << "\n";
//     }
//     fout.close();

// }

void timecalibration(TTree* cut_tree, TTree* &clonedTree){
    // Open and read .csv file
    TString dir = gSystem->UnixPathName(__FILE__);
    dir.ReplaceAll("caloanalysisPipeline.cpp", "");
    dir.ReplaceAll("/./", "/");
    // Read in time calibration data and write to ROOT file
    TFile *fin = new TFile("timecalibration.root", "RECREATE");
    TTree *t = new TTree("ntuple", "Time calibration data");
    t->ReadFile("/sps/nemo/scratch/Edinburgh_student/Katherine/timecalibration.csv", "om_number/I:mean_of_fit/F:error_of_mean/F:standard_deviation/F:standard_deviation_error/F", ';');
    fin->Write();
    // Variables
    float mof; // mean of fit (calibration value for calo_tdc)
    int omnum; // om number (from calibration)
    // Map time calibration values to om_number
    t->SetBranchAddress("mean_of_fit", &mof);
    t->SetBranchAddress("om_number", &omnum);
    std::map<int, float> calibration_map;
    Long64_t m = t->GetEntries();
    for (Long64_t j = 0; j < m; j++){
        t->GetEntry(j);
        calibration_map[omnum] = mof;
    }
    fin->Close();
    // 
    std::vector<double>* calo_tdc = nullptr;
    std::vector<int>* om_number = nullptr; 
    std::vector<double>* energy = nullptr;
    // Set Branch addresses
    cut_tree->SetBranchAddress("calo_tdc", &calo_tdc);
    cut_tree->SetBranchAddress("om_number", &om_number);
    cut_tree->SetBranchAddress("energy", &energy);
    // Output file and tree
    // TFile *fout = new TFile("calibration.root", "RECREATE");
    clonedTree = cut_tree->CloneTree(0);
    std::vector<double> cal_calo_tdc; // calibrated tdc values
    clonedTree->Branch("cal_calo_tdc", &cal_calo_tdc);
    // Event loop
    Long64_t n = cut_tree->GetEntries();
    for (Long64_t i = 0; i < n; i++){
        cut_tree->GetEntry(i);
        cal_calo_tdc.clear(); // clear previous entry
        cal_calo_tdc.reserve(om_number->size());
        for (size_t k = 0; k < om_number->size(); k++){
            int current_om = om_number->at(k);
            double current_tdc = calo_tdc->at(k);
            float cal_value = calibration_map.count(current_om) ? calibration_map.at(current_om) : 0.0f;
            if (cal_value != -999){
                cal_calo_tdc.push_back(current_tdc - cal_value);
            }
            else{
                cal_calo_tdc.push_back(current_tdc);
            }
            // std::cout << current_om << ": " << cal_calo_tdc.at(k) << ", energy: " << energy->at(k) << std::endl;
            // std::cout << "calo_tdc before calibration: " << calo_tdc->at(k) << std::endl;
        }
        clonedTree->Fill();
        // clonedTree->SetBranchStatus("*", 0); // disable irrelevant branches
        // clonedTree->SetBranchStatus("OM_count", 1);
        // clonedTree->SetBranchStatus("energy", 1);
        // clonedTree->SetBranchStatus("om_number", 1);
        // clonedTree->SetBranchStatus("event_number", 1);
        // clonedTree->SetBranchStatus("cal_calo_tdc", 1);
        // clonedTree->SetBranchStatus("gamma_om_x", 1);
        // clonedTree->SetBranchStatus("gamma_om_y", 1);
        // clonedTree->SetBranchStatus("gamma_om_z", 1);
        
        
    }
    clonedTree->ResetBranchAddresses();
    
}

void caloanalysisPipeline(const char* filename){
    // Open file and retrieve tree
    TFile *f = TFile::Open(filename);
    TTree *tree = (TTree*)f->Get("Result_tree");
    // Set up variables for existing branches
    std::vector<double>* energy = nullptr;
    std::vector<int>* om_number = nullptr;
    std::vector<int>* pid = nullptr;
    std::vector<double>* first_vertex_x = nullptr;
    std::vector<double>* second_vertex_x = nullptr;
    // Set branch addresses
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("first_vertex_x", &first_vertex_x);
    tree->SetBranchAddress("second_vertex_x", &second_vertex_x);
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
        // std::cout << OM_count << std::endl;
        // Fill new tree
        // for (double j : *energy){
        //     if (j > 12.0 && OM_count > 2) continue;
        // }
        itree->Fill();
        
    }
    // Apply selection cut (E >= 10 MeV, OM_count >= 2, electron_number == 1 for track selection)
    TTree *cut_tree = itree->CopyTree("energy >= 12.0 && OM_count >= 2");
    std::cout << cut_tree->GetEntries() << std::endl;
    // Pass cut_tree to the rest of the pipeline
    TTree* cal_tree = nullptr;
    timecalibration(cut_tree, cal_tree);
    std::cout << "Time calibration completed" << std::endl;
    // omgeometry(cal_tree);
    clusterAnalysis(cal_tree);
    std::cout << "Clustering completed" << std::endl;
    // calodistPlot();
    // clusterEnergy();
    // cluster3D();
    // Memory clean-up
    // delete cut_tree;
    delete cut_tree;
    delete cal_tree;
    f->Close();
    
}
