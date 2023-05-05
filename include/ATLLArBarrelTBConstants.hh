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

//Includers from STL
//
#include <ostream>

namespace ATLLArBarrelTBConstants{

    //explicit enum class for STAC sections
    //
    enum class STACSection : int {
        Front = 0,
        Middle = 1,
        Back =2
    };
    std::ostream& operator<<(std::ostream& ostream, const STACSection& section){
        switch (section) {
            case STACSection::Front:
                ostream << "Front Section";
                break;
            case STACSection::Middle:
                ostream << "Middle Section";
                break;
            case STACSection::Back:
                ostream << "Back Section";
                break;
        }   
        return ostream;
    }
 
    //Geometry parameters for hits and sensitive detector
    //
    const G4double EtaMin = 0.;              //we consider only half barrel (z>0)
    const G4double EtaMax = 1.475;           //end of half barrel (z>0)
    const G4double EtaChange = 0.8;          //at eta=0.8 the absorber changes
                                             //the absorber is divided in
                                             //A (eta<0.8) and B (eta>0.8) sides 
    const G4double RMin = 1.4*m;             //Internal radius of whole calo
    const G4double RMinSTAC = 1.5*m;         //Internal radius of STAC
    const G4double STAClength_eta0 = 470.*mm;//STAC length at eta=0
                                             //At eta=0 the STAC is 22.3X0
                                             //front (4.3X0) + middle (16.0X0) + back (2.0X0)
    const G4double X0 = 470./22.3;           //X0 in mm
    const G4double FrontDepth = 4.3*X0;      //Front section depth (mm)
    const G4double MiddleDepth = 16.0*X0;    //Middle section depth (mm)

} //ATLLArBarrelTBConstants

#endif //ATLLArBarrelTBConstants_h 1

//**************************************************

