//**************************************************
// \file ATLLArBarrelEventAction.cc
// \brief: implementation of the 
//         ATLLArBarrelEventAction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 6 May 2023
//**************************************************

//Includers from project files
//
#include "ATLLArBarrelEventAction.hh"
#include "ATLLArBarrelTBConstants.hh"
#include "ATLLArBarrelSensDet.hh"

//Includers from Geant4
//
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4Version.hh"
#if G4VERSION_NUMBER < 1100
#include "g4root.hh"  //replaced by G4AnalysisManager.h in G4 v11 and up
#else
#include "G4AnalysisManager.hh"
#endif
#include "G4SDManager.hh"

//Constructor and de-constructor
//
ATLLArBarrelEventAction::ATLLArBarrelEventAction()
    : G4UserEventAction(){

    fFrontHitsEdepVector  = std::vector<G4double>(ATLLArBarrelTBConstants::FrontHitNo, 0.);
    fMiddleHitsEdepVector = std::vector<G4double>(ATLLArBarrelTBConstants::MiddleHitNo, 0.);
    fBackHitsEdepVector   = std::vector<G4double>(ATLLArBarrelTBConstants::BackHitNo, 0.);
}

ATLLArBarrelEventAction::~ATLLArBarrelEventAction() {
}

//GetHitsCollection private method
//
ATLLArBarrelHitsCollection* ATLLArBarrelEventAction::GetHitsCollection(G4int hcID, const G4Event* event) const {
    
    auto HC = static_cast<ATLLArBarrelHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID)); 
    if ( ! HC ) {
        G4ExceptionDescription msg;
        msg << "Cannot access hitsCollection ID " << hcID; 
        G4Exception("ATLLArBarrelEventAction::GetHitsCollection()",
        "MyCode0003", FatalException, msg);
    }

    return HC;
}

//BeginOfEventAction base virtual method
//
void ATLLArBarrelEventAction::BeginOfEventAction([[maybe_unused]] const G4Event* event){
    
    //Clean fields of EventAction
    //
    for(auto& value : fFrontHitsEdepVector){value=0.;}
    for(auto& value : fMiddleHitsEdepVector){value=0.;}
    for(auto& value : fBackHitsEdepVector){value=0.;}
}

//EndOfEventAction base virtual method
//
void ATLLArBarrelEventAction::EndOfEventAction(const G4Event* event){
    
    //Get the hits collections and corresponding vector.
    //Each HitsCollection is identified by an integer,
    //the user shoudl retrieve the integer using the SDManager.
    //If the name of the HitCollection is not ambigous, we can skip
    //the name of the SD that created that HC.
    //
    auto FrontHCID  = G4SDManager::GetSDMpointer()->GetCollectionID(ATLLArBarrelSensDet::fFrontHitsCollectionName);
    auto MiddleHCID = G4SDManager::GetSDMpointer()->GetCollectionID(ATLLArBarrelSensDet::fMiddleHitsCollectionName); 
    auto BackHCID   = G4SDManager::GetSDMpointer()->GetCollectionID(ATLLArBarrelSensDet::fBackHitsCollectionName); 
    ATLLArBarrelHitsCollection* FrontHC  = GetHitsCollection(FrontHCID, event);
    ATLLArBarrelHitsCollection* MiddleHC = GetHitsCollection(MiddleHCID, event);
    ATLLArBarrelHitsCollection* BackHC   = GetHitsCollection(BackHCID, event);
    for (std::size_t n = 0; n < ATLLArBarrelTBConstants::FrontHitNo; ++n)  fFrontHitsEdepVector[n] = (*FrontHC)[n]->GetEdep();
    for (std::size_t n = 0; n < ATLLArBarrelTBConstants::MiddleHitNo; ++n) fMiddleHitsEdepVector[n]= (*MiddleHC)[n]->GetEdep();
    for (std::size_t n = 0; n < ATLLArBarrelTBConstants::BackHitNo; ++n)   fBackHitsEdepVector[n]  = (*BackHC)[n]->GetEdep();

}

//**************************************************
