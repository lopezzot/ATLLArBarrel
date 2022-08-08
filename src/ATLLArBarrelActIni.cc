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

ATLLArBarrelActIni::ATLLArBarrelActIni( )
    : G4VUserActionInitialization(){}

ATLLArBarrelActIni::~ATLLArBarrelActIni() {}

//Define Build() and BuildForMaster() methods
//
void ATLLArBarrelActIni::BuildForMaster() const {
}

void ATLLArBarrelActIni::Build() const {
    SetUserAction( new ATLLArBarrelPriGenAct );
}

//**************************************************
