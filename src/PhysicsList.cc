//**************************************************
// \file PhysicsList.cc
// \brief: implementation of PhysicsList class
//         for G4HepEm usage
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 20 February 2023
//**************************************************

//Includers from project files
//
#include "PhysicsList.hh" 
#include "PhysListHepEm.hh"

//Includers from Geant4
//
#include "G4LossTableManager.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList() :
    G4VModularPhysicsList() {
    
    //Set default cuts
    //
    G4LossTableManager::Instance();
    SetDefaultCutValue(0.7*mm);

    RegisterPhysics( new PhysListHepEm );
}

PhysicsList::~PhysicsList(){}

//**************************************************
