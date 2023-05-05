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
    const G4double EtaMin = 0.;                        //We consider only half barrel (z>0)
    const G4double EtaMax = 1.475;                     //end of half barrel (z>0)
    const G4double EtaChange = 0.8;                    //at eta=0.8 the absorber changes
                                                       //the absorber is divided in
                                                       //A (eta<0.8) and B (eta>0.8) sides 
    const G4double RMin = 1.4*m;                       //Internal radius of whole calo
    const G4double RMinSTAC = 1.5*m;                   //Internal radius of STAC
    const G4double STAClength_eta0 = 470.*mm;          //STAC length at eta=0
                                                       //At eta=0 the STAC is 22.3X0
                                                       //front (4.3X0) + middle (16.0X0) + back (2.0X0)
    const G4double X0 = 470./22.3;                     //X0 in mm
    const G4double FrontDepth = 4.3*X0;                //Front section depth (mm)
    const G4double MiddleDepth = FrontDepth + 16.0*X0; //Middle section depth (mm)

    const G4double BackDeltaEta = 0.05;                //(rad)DeltaEta at back section
                                                       //there are 16 eta bins in back section
                                                       //up to eta=0.8 (first half calo eta<0.8)
    const G4double MiddleDeltaEta = BackDeltaEta/2.;   //(rad)Eta bins in middle section
    const G4double FrontDeltaEta  = MiddleDeltaEta/8.; //(rad)Eta bins in front section

    const G4double STACDeltaPhi     = 2.*M_PI/16;      //(rad)Delta Phi STAC
    const G4double halfSTACDeltaPhi = STACDeltaPhi/2.;
    const G4double FrontDeltaPhi  = STACDeltaPhi/4;    //(rad)Delta Phi front section
    const G4double MiddleDeltaPhi = STACDeltaPhi/16;   //(rad)Delta Phi middle section
    const G4double BackDeltaPhi   = MiddleDeltaPhi;    //(rad)Delta Phi back section

    const G4int FrontHitNo  = 1024;                    //Number of hit to be allocated
    const G4int MiddleHitNo = 512;                     //for front, middle and back
    const G4int BackHitNo   = 256;                     //section (fist half calo eta<0.8)
    
    const G4int BackEtasPerRow = 16;                   //Eta cells per row back section
    const G4int MiddleEtasPerRow = BackEtasPerRow*2;   //Eta cells per row middle section
    const G4int FrontEtasPerRow = MiddleEtasPerRow*8;  //Eta cells per row front section

} //ATLLArBarrelTBConstants

#endif //ATLLArBarrelTBConstants_h 1

//**************************************************

