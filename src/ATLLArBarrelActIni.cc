//**************************************************
// \file ATLLArBarrelActIni.cc
// \brief: implementation of ATLLArBarrelActIni class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 8 August 2022
//**************************************************

//Includers from project files
//
#include "ATLLArBarrelActIni.hh"
#include "ATLLArBarrelPriGenAct.hh"
#include "ATLLArBarrelRunAction.hh"
#include "ATLLArBarrelEventAction.hh"
#include "ATLLArBarrelTrackingAction.hh"

ATLLArBarrelActIni::ATLLArBarrelActIni( )
    : G4VUserActionInitialization(){
    
    fEventAction = new ATLLArBarrelEventAction();
}

ATLLArBarrelActIni::~ATLLArBarrelActIni() {}

//Define Build() and BuildForMaster() methods
//
void ATLLArBarrelActIni::BuildForMaster() const {

    //From the base class comments:"The user should not use
    // this method to instantiate user action classes except for user
    // run action."
    //
    SetUserAction(new ATLLArBarrelRunAction(fEventAction));
}

void ATLLArBarrelActIni::Build() const {
    SetUserAction( new ATLLArBarrelPriGenAct );
    SetUserAction( fEventAction );
    SetUserAction( new ATLLArBarrelRunAction(fEventAction));
    SetUserAction( new ATLLArBarrelTrackingAction() );
}

//**************************************************
