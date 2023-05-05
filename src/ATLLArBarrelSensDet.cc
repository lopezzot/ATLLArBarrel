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

ATLLArBarrelSensDet::ATLLArBarrelSensDet(const G4String& name, const G4String& hitsCollectionName)
    : G4VSensitiveDetector(name) {
    
    collectionName.insert(hitsCollectionName);
}

void ATLLArBarrelSensDet::Initialize(G4HCofThisEvent* hitCollection) {}

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


    return true;
}

void ATLLArBarrelSensDet::EndOfEvent( G4HCofThisEvent* ){}

//**************************************************
