#include <iostream>
#include <algorithm>
#include <vector>

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

    //Loop over events
    //
    for(std::size_t evt; evt<tree->GetEntries(); evt++){
        tree->GetEntry(evt);

        //Find index of vector max (middle layer)
        //
        auto FrontMax = std::distance(Front->begin(),(std::max_element(Front->begin(),Front->end())));
        std::cout<<FrontMax<<std::endl;
    }

}
