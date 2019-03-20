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

    void SetFirstLayerX(G4double xx) { fFirstLayerX = xx; }
    void SetFirstLayerZ(G4double zz) { fFirstLayerZ = zz; }
    void SetFirstLayerThickness(G4double thickness) { fFirstLayerThickness = thickness; }
    void SetFirstLayerDistFromBeam(G4double dist) { fFirstLayerDistFromBeam = dist; }
    
    void SetSecondLayerX(G4double xx) { fSecondLayerX = xx; }
    void SetSecondLayerZ(G4double zz) { fSecondLayerZ = zz; }
    void SetSecondLayerThickness(G4double thickness) { fSecondLayerThickness = thickness; }
    void SetSecondLayerDistFromBeam(G4double dist) { fSecondLayerDistFromBeam = dist; }
    
    void SetThirdLayerX(G4double xx) { fThirdLayerX = xx; }
    void SetThirdLayerZ(G4double zz) { fThirdLayerZ = zz; }
    void SetThirdLayerThickness(G4double thickness) { fThirdLayerThickness = thickness; }
    void SetThirdLayerDistFromBeam(G4double dist) { fThirdLayerDistFromBeam = dist; }

private:
    
    G4int BuildTISTAR();
    
    // Assembly volumes
    G4AssemblyVolume* fAssemblyTISTAR;            

    // Logical volumes
    G4LogicalVolume * fFirstLayerLV;
    G4LogicalVolume * fSecondLayerLV;
    G4LogicalVolume * fThirdLayerLV;

    // Materials 
    G4String fSiliconMaterialName;

    // Dimensions
    G4double fFirstLayerX;
    G4double fFirstLayerZ;
    G4double fFirstLayerThickness;
    G4double fFirstLayerDistFromBeam;
       
    G4double fSecondLayerX;
    G4double fSecondLayerZ;
    G4double fSecondLayerThickness;
    G4double fSecondLayerDistFromBeam;
    
    G4double fThirdLayerX;
    G4double fThirdLayerZ;
    G4double fThirdLayerThickness;
    G4double fThirdLayerDistFromBeam;

};

#endif

