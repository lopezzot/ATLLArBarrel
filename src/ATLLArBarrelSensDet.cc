//**************************************************
// \file ATLLArBarrelSensDet.cc
// \brief: implementation of ATLLArBarrelSensDet
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 21 february 2023
//**************************************************

//Includers from project files
//
#include "ATLLArBarrelSensDet.hh"
#include "ATLLArBarrelTBConstants.hh"

//Includers from Geant4
//
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//Includers from STL
//
#include <cmath>

ATLLArBarrelSensDet::ATLLArBarrelSensDet(const G4String& name)
    : G4VSensitiveDetector(name),
    fFrontHitsCollection(nullptr),
    fMiddleHitsCollection(nullptr),
    fBackHitsCollection(nullptr){
    
    //Insert in vector of collection names the 3 collections' names
    collectionName.insert(fFrontHitsCollectionName);
    collectionName.insert(fMiddleHitsCollectionName);
    collectionName.insert(fBackHitsCollectionName);

}

void ATLLArBarrelSensDet::Initialize(G4HCofThisEvent* HCE){

    /*This method is called at the beginning of each event.
      The user must allocate with "new" the hits collections of this SD
      and link them to the SD field pointers.
      The user should as well add them to the hit collection of this event
      in such a way hit collections are automatically cleaned at the end
      of each event and can be accessed in the EndOfEvent method with
      evet->GetHCofThisEvent().*/

    //Hits Collections must have a unique pair of string (SDname + hitscollectionname)
    //
    fFrontHitsCollection  = new ATLLArBarrelHitsCollection(GetName(), fFrontHitsCollectionName); 
    fMiddleHitsCollection = new ATLLArBarrelHitsCollection(GetName(), fMiddleHitsCollectionName); 
    fBackHitsCollection   = new ATLLArBarrelHitsCollection(GetName(), fBackHitsCollectionName); 
    
    //Hits Collections (among all SD of the simulations) are also identified by an integer.
    //It can be accessed via the SDManager anytime for both storing and retrieving the hits collection.
    //If the hitcollection name is unambigous, it is enough to retrieve the ID (no need for SD name).
    //
    static G4int FrontHCID = -1, MiddleHCID = -1, BackHCID = -1;
    if(FrontHCID<0){ //prevent doing it every event -> HCIDs are identical for every event
        FrontHCID  = G4SDManager::GetSDMpointer()->GetCollectionID(fFrontHitsCollectionName); 
        MiddleHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fMiddleHitsCollectionName); 
        BackHCID   = G4SDManager::GetSDMpointer()->GetCollectionID(fBackHitsCollectionName); 
    }
    HCE->AddHitsCollection( FrontHCID, fFrontHitsCollection ); 
    HCE->AddHitsCollection( MiddleHCID, fMiddleHitsCollection ); 
    HCE->AddHitsCollection( BackHCID, fBackHitsCollection ); 

}

G4bool ATLLArBarrelSensDet::ProcessHits(G4Step* aStep, G4TouchableHistory* th){

    //Print out some info step-by-step in sensitive elements
    //
    //G4cout<<"Track #: "<< aStep->GetTrack()->GetTrackID()<< " " <<
    //        "Step #: " << aStep->GetTrack()->GetCurrentStepNumber()<< " "<<
    //        "Volume: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()<< " " << G4endl;
    //G4cout<<"x: "<< aStep->GetPreStepPoint()->GetPosition().x() <<
    //        "y: "<< aStep->GetPreStepPoint()->GetPosition().y() <<
    //        "z: "<< aStep->GetPreStepPoint()->GetPosition().z() << G4endl;
    //G4cout<<"Particle "<< aStep->GetTrack()->GetParticleDefinition()->GetParticleName()<< " " <<
    //        "Dep(MeV) "<< aStep->GetTotalEnergyDeposit() << " " <<
    //        "Mat "     << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " " << G4endl; 
    
    //Get hit Eta, Theta, Phi and Magnitude values
    //
    auto Eta = aStep->GetPreStepPoint()->GetPosition().eta();
    auto Theta = aStep->GetPreStepPoint()->GetPosition().theta(); //rad
    auto Phi = aStep->GetPreStepPoint()->GetPosition().phi();     //rad
    auto Edep = aStep->GetTotalEnergyDeposit();
    auto R = aStep->GetPreStepPoint()->GetPosition().mag();       //mm
    //If Edep==0 or Eta<0. or Eta>0.8 (i.e. only A secion allowed) do not process hit
    //
    //if(/*Edep == 0.*/false || Eta<ATLLArBarrelTBConstants::EtaMin || Eta>ATLLArBarrelTBConstants::EtaChange) return false;
    
    //Calculate depth of hit inside STAC
    //
    const G4double Depth = R - (ATLLArBarrelTBConstants::RMinSTAC/std::sin(Theta));
    ATLLArBarrelTBConstants::STACSection Section;
    G4double DeltaEta = 0.;
    G4double DeltaPhi = 0.;
    if (Depth<ATLLArBarrelTBConstants::FrontDepth){
        Section  = ATLLArBarrelTBConstants::STACSection::Front;
        DeltaEta = ATLLArBarrelTBConstants::FrontDeltaEta;
        DeltaPhi = ATLLArBarrelTBConstants::FrontDeltaPhi;
    } 
    else if (Depth<ATLLArBarrelTBConstants::MiddleDepth){
        Section  = ATLLArBarrelTBConstants::STACSection::Middle;
        DeltaEta = ATLLArBarrelTBConstants::MiddleDeltaEta;
        DeltaPhi = ATLLArBarrelTBConstants::MiddleDeltaPhi;
    }
    else{
        Section  = ATLLArBarrelTBConstants::STACSection::Back;
        DeltaEta = ATLLArBarrelTBConstants::BackDeltaEta;
        DeltaPhi = ATLLArBarrelTBConstants::BackDeltaPhi;
    }

    //Calculate Eta and Phi index and Hitmap key
    //max EtaIdx value = 256, max PhiIdx value = 16
    G4int EtaIdx = std::floor(Eta/DeltaEta);
    G4int PhiIdx = std::floor((Phi+ATLLArBarrelTBConstants::halfSTACDeltaPhi)/DeltaPhi);
    G4int HitKey = 100*EtaIdx+PhiIdx;


    G4cout<<Section<<" "<<EtaIdx<<" "<<Phi<<" "<<PhiIdx<<" "<<HitKey<<G4endl;

    return true;
}

void ATLLArBarrelSensDet::EndOfEvent( G4HCofThisEvent* ){}

//**************************************************
