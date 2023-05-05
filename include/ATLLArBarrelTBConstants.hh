//**************************************************
// \file ATLLArBarrelTBConstants.hh
// \brief: definition of ATLLArBarrelTBConstants 
//         namespace
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 5 May 2023
//**************************************************

#ifndef ATLLArBarrelTBConstants_h
#define ATLLArBarrelTBConstants_h 1

//Includers from Geant4
//
#include "G4SystemOfUnits.hh"

namespace ATLLArBarrelTBConstants{

    //Geometry parameters for hits and sensitive detector
    //
    const G4double EtaMin = 0.;     //we consider only half barrel (z>0)
    const G4double EtaMax = 1.475;  //end of half barrel (z>0)
    const G4double EtaChange = 0.8; //at eta=0.8 the absorber changes
                                    //the absorber is divided in
                                    //A (eta<0.8) and B (eta>0.8) sides 
    const G4double RMin = 1.4*m;    //Internal radius

}

#endif //ATLLArBarrelTBConstants_h 1

//**************************************************

