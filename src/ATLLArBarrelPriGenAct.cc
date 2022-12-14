//**************************************************
// \file ATLLArBarrelPriGenAct.cc
// \brief: implementation of ATLLArBarrelPriGenAct
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 8 August 2022
//**************************************************

//Includers from project files
//
#include "ATLLArBarrelPriGenAct.hh"

//Includers from Geant4
//
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

ATLLArBarrelPriGenAct::ATLLArBarrelPriGenAct()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun( nullptr ) {
    
      fParticleGun = new G4ParticleGun(1); //set primary particle(s) to 1

      //default particle gun parameters (can be changed via UI)
      //
      auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle( "e-" );
      fParticleGun->SetParticleDefinition( particleDefinition );
      fParticleGun->SetParticleMomentumDirection( G4ThreeVector( 0.,0.,1. ) );
      fParticleGun->SetParticleEnergy( 1.*GeV );
      fParticleGun->SetParticlePosition( G4ThreeVector( 0.,0.,0. ) );
}

ATLLArBarrelPriGenAct::~ATLLArBarrelPriGenAct() {
    delete fParticleGun;
}

void ATLLArBarrelPriGenAct::GeneratePrimaries( G4Event* event){
    fParticleGun->GeneratePrimaryVertex( event );
}

//**************************************************

