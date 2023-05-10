#include <iostream>
#include <algorithm>
#include <vector>

const int MiddleEtasPerRow = 16*2;

void ATLLArBarrelAnalysis(){
    
    //Open file and get tree
    //
    std::string filename = "/Users/lorenzopezzotti/Desktop/ATLBarrel/test/buildg41110/ATLLArBarrelOut_Run0.root";
    TFile* file = TFile::Open( filename.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(file->Get("ATLLArBarrelout"));

    //Set branch adresses
    //
    std::vector<double>* Front = new std::vector<double>(); 
    std::vector<double>* Middle = new std::vector<double>(); 
    std::vector<double>* Back = new std::vector<double>(); 
    tree->SetBranchAddress( "FrontHitsEdep", &Front );
    tree->SetBranchAddress( "MiddleHitsEdep",&Middle );
    tree->SetBranchAddress( "BackHitsEdep",  &Back );

    //Allocate histograms
    //
    TH1F* RPhiH1 = new TH1F("RPhi_H1","RPhi",300,0.8,1.1);

    //Loop over events
    //
    for(std::size_t evt; evt<tree->GetEntries(); evt++){
        tree->GetEntry(evt);
 
        double RPhi = 0.; //RPhi variable E3x3/E3x7 (middle layer)
 
        //Find index of vector max (middle layer)
        //
        auto MiddleMax = std::distance(Middle->begin(),(std::max_element(Middle->begin(),Middle->end())));
        //std::cout<<MiddleMax<<std::endl;
        if(MiddleMax!=271) continue; //exclude events where max is not centered in hit 271

        int MMRawM1 = MiddleMax - 1*MiddleEtasPerRow;
        int MMRawM2 = MiddleMax - 2*MiddleEtasPerRow;
        int MMRawM3 = MiddleMax - 3*MiddleEtasPerRow;
        int MMRawP1 = MiddleMax + 1*MiddleEtasPerRow;
        int MMRawP2 = MiddleMax + 2*MiddleEtasPerRow;
        int MMRawP3 = MiddleMax + 3*MiddleEtasPerRow;

        double E3x3 = Middle->at(MiddleMax)+Middle->at(MiddleMax-1)+Middle->at(MiddleMax+1)+
                      Middle->at(MMRawM1)+Middle->at(MMRawM1-1)+Middle->at(MMRawM1+1)+
                      Middle->at(MMRawP1)+Middle->at(MMRawP1-1)+Middle->at(MMRawP1+1);
                      
        double E3x7 = Middle->at(MiddleMax)+Middle->at(MiddleMax-1)+Middle->at(MiddleMax+1)+
                      Middle->at(MMRawM1)+Middle->at(MMRawM1-1)+Middle->at(MMRawM1+1)+
                      Middle->at(MMRawM2)+Middle->at(MMRawM2-1)+Middle->at(MMRawM2+1)+
                      Middle->at(MMRawM3)+Middle->at(MMRawM3-1)+Middle->at(MMRawM3+1)+
                      Middle->at(MMRawP1)+Middle->at(MMRawP1-1)+Middle->at(MMRawP1+1)+
                      Middle->at(MMRawP2)+Middle->at(MMRawP2-1)+Middle->at(MMRawP2+1)+
                      Middle->at(MMRawP3)+Middle->at(MMRawP3-1)+Middle->at(MMRawP3+1);

        RPhi = E3x3/E3x7;
        RPhiH1->Fill(RPhi);

    } //end for loop over events
    
    //Open output file and write
    //
    TFile* OutFile = TFile::Open("AnalysisBarrel.root", "RECREATE");
    RPhiH1->Write();
    OutFile->Close();

}
