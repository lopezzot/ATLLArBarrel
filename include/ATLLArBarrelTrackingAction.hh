//**************************************************
// \file ATLLArBarrelTrackingAction.hh
// \brief: definition of ATLLArBarrelTrackingAction
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 8 May 2023
//**************************************************

#ifndef ATLLArBarrelTrackingAction_h
#define ATLLArBarrelTrackingAction_h 1

//Includers from Geant4
//
#include "G4UserTrackingAction.hh"

class ATLLArBarrelTrackingAction : public G4UserTrackingAction {

    public:
        ATLLArBarrelTrackingAction();
        ~ATLLArBarrelTrackingAction();

    public:
        //Methods from base class
        //
        virtual void PreUserTrackingAction(const G4Track*);
        virtual void PostUserTrackingAction(const G4Track*);

};

#endif //ATLLArBarrelTrackingAction_h 

//**************************************************
