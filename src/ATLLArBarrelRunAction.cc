//**************************************************
// \file ATLLArBarrelRunAction.cc
// \brief: implementation of ATLLArBarrelRunAction
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 31 October 2022
//**************************************************

//Includers from project files
//
#include "ATLLArBarrelRunAction.hh"

ATLLArBarrelRunAction::ATLLArBarrelRunAction(){}

ATLLArBarrelRunAction::~ATLLArBarrelRunAction(){}

void ATLLArBarrelRunAction::BeginOfRunAction(const G4Run*){
    fTimer.Start();
}

void ATLLArBarrelRunAction::EndOfRunAction(const G4Run* aRun){
    fTimer.Stop();
    
    G4int events = aRun->GetNumberOfEvent();
    G4cout << " ============================================================================== " << G4endl;
    G4cout << "  Run terminated, " << events << " events transported" << G4endl;
    G4cout << "  Time: " << fTimer << G4endl;
    G4cout << " ============================================================================== " << G4endl;
}

//**************************************************
