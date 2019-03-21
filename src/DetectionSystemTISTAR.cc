#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemTISTAR.hh"

#include "G4SystemOfUnits.hh"

#include <string>

DetectionSystemTISTAR::DetectionSystemTISTAR() :
    fAssemblyTISTAR(NULL),
    fFirstLayerLV(NULL),
    fSecondLayerLV(NULL),
    fThirdLayerLV(NULL)
{
    // Default parameters taken from image of Leila's G4 setup
    fFirstLayerX = 40.0*mm;
    fFirstLayerZ = 50.0*mm;
    fFirstLayerThickness = 20.0*um;
    fFirstLayerDistFromBeam = 10.0*mm;    
    fFirstLayerGapZ = 4.0*mm;

    fSecondLayerX = 100.0*mm;
    fSecondLayerZ = 100.0*mm;
    fSecondLayerThickness = 150.0*um;
    fSecondLayerDistFromBeam = 30.0*mm;    
    
    fThirdLayerX = 100.0*mm;
    fThirdLayerZ = 100.0*mm;
    fThirdLayerThickness = 2000.0*um;
    fThirdLayerDistFromBeam = 33.0*mm;    
 
    fSiliconMaterialName = "Silicon";
   
}

DetectionSystemTISTAR::~DetectionSystemTISTAR() 
{
    delete fFirstLayerLV;
    delete fSecondLayerLV;
    delete fThirdLayerLV;

    delete fAssemblyTISTAR;
}

G4int DetectionSystemTISTAR::Build() 
{
    fAssemblyTISTAR = new G4AssemblyVolume(); 

    BuildTISTAR();

    return 1;
}

G4int DetectionSystemTISTAR::PlaceDetector(G4LogicalVolume* expHallLog) 
{
    G4RotationMatrix * rotate = new G4RotationMatrix;
    G4ThreeVector move = G4ThreeVector(0., 0., 0.);

    fAssemblyTISTAR->MakeImprint(expHallLog, move, rotate);

    return 1;
}

G4int DetectionSystemTISTAR::BuildTISTAR() 
{
    G4ThreeVector move, direction;
    G4RotationMatrix* rotate = new G4RotationMatrix();

    // Get the materials (should be built in the DetectorConstructionSuppressed)
    G4Material* siliconMaterial = G4Material::GetMaterial(fSiliconMaterialName);
    if( !siliconMaterial ) {
        G4cout << " ----> Material " << fSiliconMaterialName << " not found, cannot build! " << G4endl;
        return 0;
    }
    
    // Set up colours and other vis. attributes
    G4VisAttributes * siliconVisAtt = new G4VisAttributes(G4Colour::Cyan());
    siliconVisAtt->SetVisibility(true);

    // Build the first layer
    G4Box * firstLayerPV = new G4Box("firstLayerPV",fFirstLayerX/2.,fFirstLayerThickness/2.,fFirstLayerZ/2.);
    if(fFirstLayerLV == NULL) {
        fFirstLayerLV = new G4LogicalVolume(firstLayerPV,siliconMaterial,"firstLayerLV",0,0,0);
        fFirstLayerLV->SetVisAttributes(siliconVisAtt);
    }
    // Add first layer forward left strip
    move = G4ThreeVector(0.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,+fFirstLayerZ/2.+fFirstLayerGapZ/2.);
    fAssemblyTISTAR->AddPlacedVolume(fFirstLayerLV,move,rotate);
    // Add first layer backward left strip
    move = G4ThreeVector(0.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,-fFirstLayerZ/2.-fFirstLayerGapZ/2.);
    fAssemblyTISTAR->AddPlacedVolume(fFirstLayerLV,move,rotate);
    // Add first layer forward right strip
    move = G4ThreeVector(0.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,+fFirstLayerZ/2.+fFirstLayerGapZ/2.);
    fAssemblyTISTAR->AddPlacedVolume(fFirstLayerLV,move,rotate);
    // Add first layer forward right strip
    move = G4ThreeVector(0.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,-fFirstLayerZ/2.-fFirstLayerGapZ/2.);
    fAssemblyTISTAR->AddPlacedVolume(fFirstLayerLV,move,rotate);
    
    // Build the second layer
    G4Box * secondLayerPV = new G4Box("secondLayerPV",fSecondLayerX/2.,fSecondLayerThickness/2.,fSecondLayerZ/2.);
    if(fSecondLayerLV == NULL) {
        fSecondLayerLV = new G4LogicalVolume(secondLayerPV,siliconMaterial,"secondLayerLV",0,0,0);
        fSecondLayerLV->SetVisAttributes(siliconVisAtt);
    }
    // Add second layer left strip
    move = G4ThreeVector(0.,+fSecondLayerThickness/2.+fSecondLayerDistFromBeam,0.);
    fAssemblyTISTAR->AddPlacedVolume(fSecondLayerLV,move,rotate);
    // Add second layer right strip
    move = G4ThreeVector(0.,-fSecondLayerThickness/2.-fSecondLayerDistFromBeam,0.);
    fAssemblyTISTAR->AddPlacedVolume(fSecondLayerLV,move,rotate);

    // Build the third layer
    G4Box * thirdLayerPV = new G4Box("thirdLayerPV",fThirdLayerX/2.,fThirdLayerThickness/2.,fThirdLayerZ/2.);
    if(fThirdLayerLV == NULL) {
        fThirdLayerLV = new G4LogicalVolume(thirdLayerPV,siliconMaterial,"thirdLayerLV",0,0,0);
        fThirdLayerLV->SetVisAttributes(siliconVisAtt);
    }
    // Add third left strip
    move = G4ThreeVector(0.,+fThirdLayerThickness/2.+fThirdLayerDistFromBeam,0.);
    fAssemblyTISTAR->AddPlacedVolume(fThirdLayerLV,move,rotate);
    // Add third right strip
    move = G4ThreeVector(0.,-fThirdLayerThickness/2.-fThirdLayerDistFromBeam,0.);
    fAssemblyTISTAR->AddPlacedVolume(fThirdLayerLV,move,rotate);

    return 1;
}

