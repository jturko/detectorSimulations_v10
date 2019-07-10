//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DETECTIONSYSTEMTISTAR_HH
#define DETECTIONSYSTEMTISTAR_HH

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

class DetectionSystemTISTAR
{
public:
    DetectionSystemTISTAR();
    ~DetectionSystemTISTAR();

    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* expHallLog);
    G4int PlaceDetector(G4ThreeVector move, G4ThreeVector rotate, G4LogicalVolume* expHallLog);

    void SetSiDimensions(G4ThreeVector dim) { fSiDimensions = dim; fSiDimensionsSet = true; }
    void SetPCBDimensions(G4ThreeVector dim) { fPCBDimensions = dim; fPCBDimensionsSet = true; }
    void SetSiOffsetInPCB(G4ThreeVector offset) { fSiOffsetInPCB = offset; } // this is the offset of the Si layer within the PCB, probably should be re-named at some point...
    void SetPositionOffset(G4ThreeVector offset) { fPositionOffset = offset; }

    G4int BuildLayer();

    G4int Add2StripLayer(G4double dist_from_beam, G4bool si_centered, G4LogicalVolume* expHallLog);
    G4int Add4StripLayer(G4double dist_from_beam, G4double gap_z, G4LogicalVolume* expHallLog);
    
    G4int AddGasTarget(G4LogicalVolume* expHallLog);

    void SetDetectorNumber(G4int detNum) { fDetectorNumber = detNum; }

    // Vacuum chamber
    void SetVacuumChamberShape(G4String shape);
    void SetVacuumChamberMaterialName(G4String material) { fVacuumChamberMaterialName = material; }
    void SetVacuumChamberBoxDimensions(G4ThreeVector dims) { fVacuumChamberBoxDimensions = dims; }
    void SetVacuumChamberCylinderRadius(G4double radius) { fVacuumChamberCylinderRadius = radius; }
    void SetVacuumChamberCylinderZ(G4double zz) { fVacuumChamberCylinderZ = zz; }
    
    void SetVacuumChamberExteriorMaterial(G4String material) { fVacuumChamberExteriorMaterialName = material; }
    void SetVacuumChamberExteriorThickness(G4double thickness) { fVacuumChamberExteriorThickness = thickness; }    

    G4int AddVacuumChamber(G4LogicalVolume* expHallLog, G4LogicalVolume *& vaccumChamberLog);

private:
    G4int fDetectorNumber;
    
    // Assembly volumes
    G4AssemblyVolume* fAssemblyLayer;

    // Logical volumes
    G4LogicalVolume* fLogicalSiLayer;
    G4LogicalVolume* fLogicalPCBLayer;

    // Materials 
    G4String fSiliconMaterialName;
    G4String fPCBMaterialName;

    // Dimensions
    G4ThreeVector fSiDimensions;
    G4ThreeVector fPCBDimensions;
    G4ThreeVector fSiOffsetInPCB;
    G4ThreeVector fPositionOffset;
    G4bool fSiDimensionsSet; 
    G4bool fPCBDimensionsSet; 
        
    // Gas target
    G4LogicalVolume* fLogicalGasTarget;
    G4double fGasTargetRadius;
    G4double fGasTargetLength;
    G4double fGasTargetDensity;
    G4double fGasTargetPressure;
    G4String fGasTargetMaterialName;

    G4LogicalVolume* fLogicalGasTargetMylar;
    G4double fGasTargetMylarThickness;
    G4String fGasTargetMylarMaterialName;

    G4LogicalVolume* fLogicalGasTargetBeWindow;
    G4double fGasTargetBeWindowThickness;
    G4String fGasTargetBeWindowMaterialName;

    // Vacuum chamber
    G4String fVacuumChamberShape;
    G4String fVacuumChamberMaterialName;
    G4ThreeVector fVacuumChamberBoxDimensions;
    G4double fVacuumChamberCylinderRadius;
    G4double fVacuumChamberCylinderZ;

    G4String fVacuumChamberExteriorMaterialName;
    G4double fVacuumChamberExteriorThickness; 
           
};

#endif
