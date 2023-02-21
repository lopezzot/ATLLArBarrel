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

//Includers from Geant4
//
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

ATLLArBarrelSensDet::ATLLArBarrelSensDet(const G4String& name, const G4String& hitsCollectionName)
    : G4VSensitiveDetector(name) {
    
    collectionName.insert(hitsCollectionName);
}

void ATLLArBarrelSensDet::Initialize(G4HCofThisEvent* hitCollection) {}

G4bool ATLLArBarrelSensDet::ProcessHits(G4Step* step, G4TouchableHistory* th){}

void ATLLArBarrelSensDet::EndOfEvent( G4HCofThisEvent* ){}

//**************************************************
