#include <iostream>
#include <algorithm>
#include <vector>

void ATLLArBarrelAnalysis(){
    
    //Open file and get tree
    //
    //std::string filename = "/Users/lorenzopezzotti/Desktop/ATLBarrel/test/buildg41110/ATLLArBarrelOut_Run0.root";
    //std::string filename = "/Users/lorenzopezzotti/Desktop/ATLBarrel/test/buildg41110/ATLLArBarrelOut_Run0_35e-.root";
    //std::string filename = "/Users/lorenzopezzotti/Desktop/ATLBarrel/test/buildg41110/eta0.2phi0.01225/20gamma.root";
    //std::string filename = "/Users/lorenzopezzotti/Desktop/ATLBarrel/test/buildg41110/eta0.2phi0.0245_noise/20gamma.root";
    //std::string filename = "/Users/lorenzopezzotti/Desktop/ATLBarrel/test/buildg41110/eta0.2phi0.0245_noise_0.1/20gamma.root";
    std::string filename = "/Users/lorenzopezzotti/Desktop/ATLBarrel/test/buildg41110/eta0.2phi0.0245_noise_0.01/20gamma.root";
    TFile* file = TFile::Open( filename.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(file->Get("ATLLArBarrelout"));

    const int MiddleEtasPerRow = 16*2;
 
    //Set branch adresses
    //
    std::vector<double>* Front = nullptr; 
    std::vector<double>* Middle = nullptr; 
    std::vector<double>* Back = nullptr; 
    tree->SetBranchAddress( "FrontHitsEdep", &Front );
    tree->SetBranchAddress( "MiddleHitsEdep",&Middle );
    tree->SetBranchAddress( "BackHitsEdep",  &Back );

    //Allocate histograms
    //
    TH1F* RPhiH1 = new TH1F("RPhi_H1","RPhi",20,0.7,1.1);
    TH1F* REtaH1 = new TH1F("REta_H1","REta",20,0.7,1.1);

    //Loop over events
    //
    for(int evt=0; evt<tree->GetEntries(); evt++){
        tree->GetEntry(evt);
 
        double RPhi = 0.; //RPhi variable E3x3/E3x7 (middle layer)
        double REta = 0.; //REta variable E3x3/E7x7 (middle layer)
 
        //Find index of vector max (middle layer)
        //
        auto MiddleMax = std::distance(Middle->begin(),(std::max_element(Middle->begin(),Middle->end())));
        //std::cout<<MiddleMax<<std::endl;
        //if(MiddleMax!=271) continue; //exclude events where max is not centered in hit 271

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

        double E7x7 = Middle->at(MiddleMax)+Middle->at(MiddleMax-1)+Middle->at(MiddleMax+1)+Middle->at(MiddleMax-2)+Middle->at(MiddleMax+2)+Middle->at(MiddleMax-3)+Middle->at(MiddleMax+3)+
                      Middle->at(MMRawM1)+Middle->at(MMRawM1-1)+Middle->at(MMRawM1+1)+Middle->at(MMRawM1-2)+Middle->at(MMRawM1+2)+Middle->at(MMRawM1-3)+Middle->at(MMRawM1+3)+
                      Middle->at(MMRawM2)+Middle->at(MMRawM2-1)+Middle->at(MMRawM2+1)+Middle->at(MMRawM2-2)+Middle->at(MMRawM2+2)+Middle->at(MMRawM2-3)+Middle->at(MMRawM2+3)+
                      Middle->at(MMRawM3)+Middle->at(MMRawM3-1)+Middle->at(MMRawM3+1)+Middle->at(MMRawM3-2)+Middle->at(MMRawM3+2)+Middle->at(MMRawM3-3)+Middle->at(MMRawM3+3)+
                      Middle->at(MMRawP1)+Middle->at(MMRawP1-1)+Middle->at(MMRawP1+1)+Middle->at(MMRawP1-2)+Middle->at(MMRawP1+2)+Middle->at(MMRawP1-3)+Middle->at(MMRawP1+3)+
                      Middle->at(MMRawP2)+Middle->at(MMRawP2-1)+Middle->at(MMRawP2+1)+Middle->at(MMRawP2-2)+Middle->at(MMRawP2+2)+Middle->at(MMRawP2-3)+Middle->at(MMRawP2+3)+
                      Middle->at(MMRawP3)+Middle->at(MMRawP3-1)+Middle->at(MMRawP3+1)+Middle->at(MMRawP3-2)+Middle->at(MMRawP3+2)+Middle->at(MMRawP3-3)+Middle->at(MMRawP3+3);

        RPhi = E3x3/E3x7;
        REta = E3x7/E7x7;
        RPhiH1->Fill(RPhi);
        REtaH1->Fill(REta);

    } //end for loop over events
 
    //Open output file and write
    //
    TFile* OutFile = TFile::Open("OutputAnalysisBarrel.root", "RECREATE");
    RPhiH1->Write();
    REtaH1->Write();
    OutFile->Close();
}
