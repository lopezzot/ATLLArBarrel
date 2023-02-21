//**************************************************
// \file ATLLArBarrelSensDet.hh
// \brief: definition of ATLLArBarrelSensDet class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 21 february 2023
//**************************************************

#ifndef ATLLArBarrelSensDet_h
#define ATLLArBarrelSensDet_h 1

//Includers from Geant4
//
#include "G4VSensitiveDetector.hh"

//Forward declarations from Geant4
//
class G4Step;
class G4HCofThisEvent;

//Includers from STL
//
#include <vector>

class ATLLArBarrelSensDet : public G4VSensitiveDetector {

    public:
        ATLLArBarrelSensDet(const G4String& name, const G4String& hitsCollectionName);
        ~ATLLArBarrelSensDet() override = default;

        //methods from base class
        //
        void   Initialize(G4HCofThisEvent* hitCollection) override;
        G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
        void   EndOfEvent(G4HCofThisEvent* hitCollection) override;

    private:
};

#endif

//**************************************************
