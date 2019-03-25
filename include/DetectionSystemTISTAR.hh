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
    DetectionSystemTISTAR(G4int nlayers=1);
    ~DetectionSystemTISTAR();

    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* expHallLog);
    G4int PlaceDetector(G4int layer, G4ThreeVector move, G4ThreeVector rotate, G4LogicalVolume* expHallLog);

    void SetDimensionsSiLayers(G4int i, G4ThreeVector dim) { fDimensionsSiLayers.at(i) = dim; fDimensionsSet.at(i).at(0) = true; }
    void SetDimensionsPCBLayers(G4int i, G4ThreeVector dim) { fDimensionsPCBLayers.at(i) = dim; fDimensionsSet.at(i).at(1) = true; }
    void SetOffsetLayers(G4int i, G4ThreeVector offset) { fOffsetLayers.at(i) = offset; }

    G4int BuildLayer(G4int layer);

private:
    
    G4int BuildTISTAR();
    
    G4int BuildSiliconStrips();
    G4int BuildPCBs();

    void CalculateDimensions();

    const G4int fNLayers;
    
    // Assembly volumes
    std::vector<G4AssemblyVolume*> fAssemblyLayers;
    
    // Logical volumes
    std::vector<G4LogicalVolume*> fLogicalSiLayers;
    std::vector<G4LogicalVolume*> fLogicalPCBLayers;
    
    // Materials 
    G4String fSiliconMaterialName;
    G4String fPCBMaterialName;

    // Dimensions
    std::vector<G4ThreeVector> fDimensionsSiLayers;
    std::vector<G4ThreeVector> fDimensionsPCBLayers;
    std::vector<G4ThreeVector> fOffsetLayers;
    std::vector<std::vector<G4bool>> fDimensionsSet; 

};

#endif

