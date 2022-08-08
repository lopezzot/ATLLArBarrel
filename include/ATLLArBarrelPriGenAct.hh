//**************************************************
// \file ATLLArBarrelPriGenAct.hh
// \brief: definition of ATLLArBarrelPriGenAct class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 8 August 2022
//**************************************************

#ifndef ATLLArBarrel_h
#define ATLLArBarrel_h 1

//Includers from Geant4
//
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class ATLLArBarrelPriGenAct : public G4VUserPrimaryGeneratorAction {

    public:
        ATLLArBarrelPriGenAct();
        virtual ~ATLLArBarrelPriGenAct();

        virtual void GeneratePrimaries( G4Event* event );

    private:
        G4ParticleGun* fParticleGun;
};

#endif //ATLLArBarrel_h 1

//**************************************************

