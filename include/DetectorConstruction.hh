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
/// \file analysis/shared/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.hh 77256 2013-11-22 10:10:23Z gcosmo $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
//class DetectionSystemGammaTracking;
class DetectionSystemGriffin;
class DetectionSystem8pi;

class ApparatusDescantStructure;
class ApparatusLayeredTarget;
class DetectionSystemDescant;

class DetectionSystemSceptar;
class DetectionSystemSpice;
class DetectionSystemPaces;
class DetectionSystemSodiumIodide;
class DetectionSystemLanthanumBromide;
class DetectionSystemBox;
class DetectionSystemAncillaryBGO;

class DetectionSystemTISTAR;

//class MagneticField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    
  
    DetectorConstruction();
    ~DetectorConstruction();
    
    void PassEfficiencyPosition( G4ThreeVector num ) {fDetEffPosition = num;}

    G4int fGriffinDetectorsMapIndex;
    G4int fGriffinDetectorsMap[16];

    void SetWorldMaterial( G4String );
    void SetWorldDimensions( G4ThreeVector );
    void SetWorldVis( G4bool );
    //    void SetWorldMagneticField( G4ThreeVector );

    //Gneeric Target
    void SetGenericTargetMaterial( G4String );
    void SetGenericTargetDimensions( G4ThreeVector );
    void SetGenericTargetPosition( G4ThreeVector );
    void SetGenericTarget( );
    
    G4double LayeredTargetLayerStart(int);
 
    //    void SetFieldBoxMaterial( G4String );
    //    void SetFieldBoxDimensions( G4ThreeVector );
    //    void SetFieldBoxPosition( G4ThreeVector );
    //    void SetFieldBoxMagneticField( G4ThreeVector );
    //    void SetFieldBox( );

    //    void SetBoxMat( G4String input )                   {boxMat = input;};
    //        void SetBoxThickness( G4double input )             {boxThickness = input;};
    //        void SetBoxInnerDimensions( G4ThreeVector input )  {boxInnerDimensions = input;};
    //        void SetBoxColour( G4ThreeVector input )           {boxColour = input;};
    //        void AddBox();

    void LayeredTargetAdd(G4String, G4double);
    
    void SetTabMagneticField(G4String, G4double, G4double);
    // Grid Functions
    void SetGridMat( G4String input )                  {fGridMat = input;};
    void SetGridSize( G4double input )                 {fGridSize = input;};
    void SetGridDimensions( G4ThreeVector input )      {fGridDimensions = input;};
    void SetGridColour( G4ThreeVector input )          {fGridColour = input;};
    void SetGridPosOffset( G4ThreeVector input )          {fGridOffset = input;};
    void AddGrid();
    void AddApparatusSpiceTargetChamber(G4String);
    void AddApparatus8piVacuumChamber();
    void AddApparatus8piVacuumChamberAuxMatShell(G4double thickness);
    void AddApparatusGriffinStructure(G4int selector);

    G4double GetWorldSizeX()           {return fWorldSizeX;};
    G4double GetWorldSizeY()           {return fWorldSizeY;};
    G4double GetWorldSizeZ()           {return fWorldSizeZ;};


    const G4VPhysicalVolume* GetphysiWorld() {return fPhysiWorld;};

    G4VPhysicalVolume* Construct();

    void UpdateGeometry();

    //    void AddDetectionSystemGammaTracking(G4int ndet);
    void AddDetectionSystemSodiumIodide(G4int ndet);
    void AddDetectionSystemLanthanumBromide(G4ThreeVector input);
    void AddDetectionSystemAncillaryBGO(G4ThreeVector input);

    void AddDetectionSystem8pi(G4int ndet);

    void AddDetectionSystemDescant(G4int ndet);
    void AddDetectionSystemDescantAuxPorts(G4ThreeVector input);
    void SetDetectionSystemDescantRotation(G4ThreeVector input);
    void SetDetectionSystemDescantColor(G4String input);
    void AddDetectionSystemDescantCart(G4ThreeVector input);
    void AddDetectionSystemDescantSpher(G4ThreeVector input, G4double unit);
    void AddApparatusDescantStructure();

    void AddDetectionSystemTestcan(G4ThreeVector input);

    void AddDetectionSystem8piDetector(G4int ndet);
    void AddDetectionSystemGriffinForward(G4int ndet);
    void AddDetectionSystemGriffinForwardDetector(G4int ndet);
    void AddDetectionSystemGriffinBack(G4int ndet);
    void AddDetectionSystemGriffinBackDetector(G4int ndet);
    //void AddDetectionSystemGriffinPositionConfig(G4ThreeVector input);
    void AddDetectionSystemGriffinHevimet(G4int input);
    void AddDetectionSystemGriffinCustom(G4int ndet);
    void AddDetectionSystemGriffinCustomDetector(G4int ndet  = 0);
    void AddDetectionSystemGriffinShieldSelect(G4int ShieldSelect );
    void AddDetectionSystemGriffinSetRadialDistance(G4double detectorDist);
    void AddDetectionSystemGriffinSetExtensionSuppLocation(G4int detectorPos);
    void AddDetectionSystemGriffinSetDeadLayer(G4ThreeVector params);

    void AddDetectionSystemSceptar(G4int ndet);
    void AddDetectionSystemPaces(G4int ndet);
    void AddDetectionSystemSpice(G4int nRings);

    void AddDetectionSystemTISTAR();

    G4double GetLanthanumBromideCrystalRadius();
    G4double GetLanthanumBromideCrystalLength();
    G4double GetLanthanumBromideR();
    G4double GetLanthanumBromideTheta(G4int i);
    G4double GetLanthanumBromidePhi(G4int i);
    G4double GetLanthanumBromideYaw(G4int i);
    G4double GetLanthanumBromidePitch(G4int i);
    G4double GetLanthanumBromideRoll(G4int i);
    G4double GetLanthanumBromideCrystalRadialPosition();

    void SetSpiceRes(G4bool);
    G4bool GetSpiceIn(){return fSetSpiceIn;};
    void UseTIGRESSPositions( G4bool input )                  {fUseTigressPositions = input;};
    
    void SetTISTARFirstLayerX(G4double value)               { fTISTARFirstLayerX = value; } 
    void SetTISTARFirstLayerZ(G4double value)               { fTISTARFirstLayerZ = value; } 
    void SetTISTARFirstLayerThickness(G4double value)       { fTISTARFirstLayerThickness = value; } 
    void SetTISTARFirstLayerDistFromBeam(G4double value)    { fTISTARFirstLayerDistFromBeam = value; } 
    void SetTISTARFirstLayerOffset(G4ThreeVector value)     { fTISTARFirstLayerOffset = value; } 
    void SetTISTARFirstLayerPCBUpperX(G4double value)       { fTISTARFirstLayerPCBUpperX = value; }
    void SetTISTARFirstLayerPCBLowerX(G4double value)       { fTISTARFirstLayerPCBLowerX = value; }
    void SetTISTARFirstLayerPCBForwardZ(G4double value)     { fTISTARFirstLayerPCBForwardZ = value; }
    void SetTISTARFirstLayerPCBBackwardZ(G4double value)    { fTISTARFirstLayerPCBBackwardZ = value; }
    void SetTISTARFirstLayerGapZ(G4double value)            { fTISTARFirstLayerGapZ = value; }

    void SetTISTARSecondLayerX(G4double value)              { fTISTARSecondLayerX = value; } 
    void SetTISTARSecondLayerZ(G4double value)              { fTISTARSecondLayerZ = value; } 
    void SetTISTARSecondLayerThickness(G4double value)      { fTISTARSecondLayerThickness = value; } 
    void SetTISTARSecondLayerDistFromBeam(G4double value)   { fTISTARSecondLayerDistFromBeam = value; } 
    void SetTISTARSecondLayerOffset(G4ThreeVector value)    { fTISTARSecondLayerOffset = value; } 
    void SetTISTARSecondLayerPCBUpperX(G4double value)       { fTISTARSecondLayerPCBUpperX = value; }
    void SetTISTARSecondLayerPCBLowerX(G4double value)       { fTISTARSecondLayerPCBLowerX = value; }
    void SetTISTARSecondLayerPCBForwardZ(G4double value)     { fTISTARSecondLayerPCBForwardZ = value; }
    void SetTISTARSecondLayerPCBBackwardZ(G4double value)    { fTISTARSecondLayerPCBBackwardZ = value; }

    void SetTISTARThirdLayerX(G4double value)               { fTISTARThirdLayerX = value; } 
    void SetTISTARThirdLayerZ(G4double value)               { fTISTARThirdLayerZ = value; } 
    void SetTISTARThirdLayerThickness(G4double value)       { fTISTARThirdLayerThickness = value; } 
    void SetTISTARThirdLayerDistFromBeam(G4double value)    { fTISTARThirdLayerDistFromBeam = value; } 
    void SetTISTARThirdLayerOffset(G4ThreeVector value)     { fTISTARThirdLayerOffset = value; } 
    void SetTISTARThirdLayerPCBUpperX(G4double value)       { fTISTARThirdLayerPCBUpperX = value; }
    void SetTISTARThirdLayerPCBLowerX(G4double value)       { fTISTARThirdLayerPCBLowerX = value; }
    void SetTISTARThirdLayerPCBForwardZ(G4double value)     { fTISTARThirdLayerPCBForwardZ = value; }
    void SetTISTARThirdLayerPCBBackwardZ(G4double value)    { fTISTARThirdLayerPCBBackwardZ = value; }

    void AddTISTARLayer();
    void SetTISTARSiDimensions(G4ThreeVector dim)   { fTISTARSiDimensions = dim; } 
    void SetTISTARPCBDimensions(G4ThreeVector dim)  { fTISTARPCBDimensions = dim; } 
    void SetTISTAROffset(G4ThreeVector offset)      { fTISTAROffset = offset; } 
    void SetTISTARRotation(G4ThreeVector rotate)    { fTISTARRotation = rotate; }
    void SetTISTARPosition(G4ThreeVector move)      { fTISTARPosition = move; }
    

private:

    //    MagneticField* worldMagField;

    G4double  fWorldSizeX;
    G4double  fWorldSizeY;
    G4double  fWorldSizeZ;
    G4bool    fWorldVis;
    G4bool    fBuiltDetectors;
    G4double  fGriffinFwdBackPosition;
    G4int     fDetectorShieldSelect ;
    G4double  fDetectorRadialDistance ;
    G4int     fExtensionSuppressorLocation ;
    G4int     fCustomDetectorNumber ;
    G4int     fCustomDetectorPosition ;
    G4int     fCustomDetectorVal ;
    G4int     fHevimetSelector ;
    G4bool    fUseTigressPositions;

    // Box
    G4String           fBoxMat;
    G4double           fBoxThickness;
    G4ThreeVector      fBoxInnerDimensions;
    G4ThreeVector      fBoxColour;

    G4Box*             fSolidWorld;    //pointer to the solid World
    G4LogicalVolume*   fLogicWorld;    //pointer to the logical World
    G4VPhysicalVolume* fPhysiWorld;    //pointer to the physical World

    // Grid
    G4String           fGridMat;
    G4double           fGridSize;
    G4ThreeVector      fGridDimensions;
    G4ThreeVector      fGridColour;
    G4ThreeVector      fGridOffset;

    void DefineSuppressedParameters();
    void DefineMaterials();

    G4double fCoords[20][5];
    G4bool        fSetGenericTargetMaterial;
    G4bool        fSetGenericTargetDimensions;
    G4bool        fSetGenericTargetPosition;
    G4String      fGenericTargetMaterial;
    G4ThreeVector fGenericTargetDimensions;
    G4ThreeVector fGenericTargetPosition;
        
    G4bool        fSetSpiceIn;
    
private: 
    G4bool        fSetFieldBoxMaterial;
    G4bool        fSetFieldBoxDimensions;
    G4bool        fSetFieldBoxPosition;
    G4bool        fSetFieldBoxMagneticField;
    G4String      fFieldBoxMaterial;
    G4ThreeVector fFieldBoxDimensions;
    G4ThreeVector fFieldBoxPosition;
    G4ThreeVector fFieldBoxMagneticField;

    G4String fMatWorldName;

    ApparatusLayeredTarget* fApparatusLayeredTarget;
    DetectorMessenger* fDetectorMessenger;
    

    G4ThreeVector fDescantRotation;
    G4String fDescantColor;
    
    G4ThreeVector fDetEffPosition;

    // TI-STAR Dimensions
    G4double fTISTARFirstLayerX;
    G4double fTISTARFirstLayerZ;
    G4double fTISTARFirstLayerThickness;
    G4double fTISTARFirstLayerDistFromBeam;
    G4ThreeVector fTISTARFirstLayerOffset;
    G4double fTISTARFirstLayerPCBUpperX;
    G4double fTISTARFirstLayerPCBLowerX;
    G4double fTISTARFirstLayerPCBForwardZ;
    G4double fTISTARFirstLayerPCBBackwardZ;
    G4double fTISTARFirstLayerGapZ;

    G4double fTISTARSecondLayerX;
    G4double fTISTARSecondLayerZ;
    G4double fTISTARSecondLayerThickness;
    G4double fTISTARSecondLayerDistFromBeam;
    G4ThreeVector fTISTARSecondLayerOffset;
    G4double fTISTARSecondLayerPCBUpperX;
    G4double fTISTARSecondLayerPCBLowerX;
    G4double fTISTARSecondLayerPCBForwardZ;
    G4double fTISTARSecondLayerPCBBackwardZ;

    G4double fTISTARThirdLayerX;
    G4double fTISTARThirdLayerZ;
    G4double fTISTARThirdLayerThickness;
    G4double fTISTARThirdLayerDistFromBeam;
    G4ThreeVector fTISTARThirdLayerOffset;
    G4double fTISTARThirdLayerPCBUpperX;
    G4double fTISTARThirdLayerPCBLowerX;
    G4double fTISTARThirdLayerPCBForwardZ;
    G4double fTISTARThirdLayerPCBBackwardZ;

    G4ThreeVector fTISTARSiDimensions;
    G4ThreeVector fTISTARPCBDimensions;
    G4ThreeVector fTISTAROffset;
    G4ThreeVector fTISTARRotation;
    G4ThreeVector fTISTARPosition;

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


