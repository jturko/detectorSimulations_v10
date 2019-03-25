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

DetectionSystemTISTAR::DetectionSystemTISTAR(G4int nlayers) :
    fNLayers(nlayers),
    fAssemblyTISTAR(NULL),
    fFirstLayerLV(NULL),
    fSecondLayerLV(NULL),
    fThirdLayerLV(NULL),
    fFirstLayerPCBUpperXLV(NULL),
    fFirstLayerPCBLowerXLV(NULL),
    fFirstLayerPCBForwardZLV(NULL),
    fFirstLayerPCBBackwardZLV(NULL),
    fSecondLayerPCBUpperXLV(NULL),
    fSecondLayerPCBLowerXLV(NULL),
    fSecondLayerPCBForwardZLV(NULL),
    fSecondLayerPCBBackwardZLV(NULL),
    fThirdLayerPCBUpperXLV(NULL),
    fThirdLayerPCBLowerXLV(NULL),
    fThirdLayerPCBForwardZLV(NULL),
    fThirdLayerPCBBackwardZLV(NULL)
{
    // Default parameters taken from image of Leila's G4 setup
    fFirstLayerX = 40.0*mm;
    fFirstLayerZ = 50.0*mm;
    fFirstLayerThickness = 20.0*um;
    fFirstLayerDistFromBeam = 10.0*mm;    
    fFirstLayerGapZ = 4.0*mm;

    fFirstLayerPCBUpperX = 0.0*mm;
    fFirstLayerPCBLowerX = 0.0*mm;
    fFirstLayerPCBForwardZ = 0.0*mm;
    fFirstLayerPCBBackwardZ = 0.0*mm;

    fSecondLayerX = 100.0*mm;
    fSecondLayerZ = 100.0*mm;
    fSecondLayerThickness = 150.0*um;
    fSecondLayerDistFromBeam = 30.0*mm;    

    fSecondLayerPCBUpperX = 0.0*mm;
    fSecondLayerPCBLowerX = 0.0*mm;
    fSecondLayerPCBForwardZ = 0.0*mm;
    fSecondLayerPCBBackwardZ = 0.0*mm;
    
    fThirdLayerX = 100.0*mm;
    fThirdLayerZ = 100.0*mm;
    fThirdLayerThickness = 2000.0*um;
    fThirdLayerDistFromBeam = 33.0*mm;    

    fThirdLayerPCBUpperX = 0.0*mm;
    fThirdLayerPCBLowerX = 0.0*mm;
    fThirdLayerPCBForwardZ = 0.0*mm;
    fThirdLayerPCBBackwardZ = 0.0*mm;
 
    fSiliconMaterialName = "Silicon";
    fPCBMaterialName = "Epoxy-Resin";
   
    fPCBThickness = 1.57*mm;

    fAssemblyLayers.resize(fNLayers); 
    fLogicalSiLayers.resize(fNLayers);
    fLogicalPCBLayers.resize(fNLayers);
    fDimensionsSiLayers.resize(fNLayers);
    fDimensionsPCBLayers.resize(fNLayers);
    fOffsetLayers.resize(fNLayers);
    for(int i=0; i<fNLayers; i++) {
        fAssemblyLayers.at(i) = NULL;
        fLogicalSiLayers.at(i) = NULL;
        fLogicalPCBLayers.at(i) = NULL;
    }

    fDimensionsSet.resize(fNLayers);
    for(int i=0; i<fNLayers; ++i) {
        fDimensionsSet.at(i).resize(2);
        fDimensionsSet.at(i).at(0) = false; // silicon dimensons
        fDimensionsSet.at(i).at(1) = false; // PCB dimensions
        fOffsetLayers.at(i) = G4ThreeVector(0.,0.,0.);
    }
}

DetectionSystemTISTAR::~DetectionSystemTISTAR() 
{
    delete fFirstLayerLV;
    delete fSecondLayerLV;
    delete fThirdLayerLV;

    delete fFirstLayerPCBUpperXLV;
    delete fFirstLayerPCBLowerXLV;
    delete fFirstLayerPCBForwardZLV;
    delete fFirstLayerPCBBackwardZLV;
    
    delete fSecondLayerPCBUpperXLV;
    delete fSecondLayerPCBLowerXLV;
    delete fSecondLayerPCBForwardZLV;
    delete fSecondLayerPCBBackwardZLV;
    
    delete fThirdLayerPCBUpperXLV;
    delete fThirdLayerPCBLowerXLV;
    delete fThirdLayerPCBForwardZLV;
    delete fThirdLayerPCBBackwardZLV;

    delete fAssemblyTISTAR;

    for(int i=0; i<int(fLogicalSiLayers.size()); ++i)   delete (fLogicalSiLayers[i]);
    for(int i=0; i<int(fLogicalPCBLayers.size()); ++i)  delete (fLogicalPCBLayers[i]);
    for(int i=0; i<int(fAssemblyLayers.size()); ++i)    delete (fAssemblyLayers[i]);
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

    //fAssemblyTISTAR->MakeImprint(expHallLog, move, rotate);
    fAssemblyLayers.at(0)->MakeImprint(expHallLog, move, rotate);

    return 1;
}

G4int DetectionSystemTISTAR::PlaceDetector(G4int layer, G4ThreeVector move, G4ThreeVector rotate, G4LogicalVolume* expHallLog)
{
    if(layer<0||layer>=fNLayers) {
        G4cout << " ---> This layer doesn't exist! " << G4endl;
        return 0;
    }

    G4RotationMatrix * rotation = new G4RotationMatrix;
    rotation->rotateX(rotate.x()*M_PI/180.);
    rotation->rotateY(rotate.y()*M_PI/180.);
    rotation->rotateZ(rotate.z()*M_PI/180.);
    fAssemblyLayers.at(layer)->MakeImprint(expHallLog, move, rotation);

    return 1;
}

G4int DetectionSystemTISTAR::BuildTISTAR() 
{

    //BuildSiliconStrips();
    //BuildPCBs();
    BuildLayer(0);

    return 1;
}

G4int DetectionSystemTISTAR::BuildSiliconStrips()
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
    move = G4ThreeVector(0.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,+fFirstLayerZ/2.+fFirstLayerGapZ/2.) + fFirstLayerOffset;
    fAssemblyTISTAR->AddPlacedVolume(fFirstLayerLV,move,rotate);
    // Add first layer backward left strip
    move = G4ThreeVector(0.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,-fFirstLayerZ/2.-fFirstLayerGapZ/2.) + fFirstLayerOffset;
    fAssemblyTISTAR->AddPlacedVolume(fFirstLayerLV,move,rotate);
    // Add first layer forward right strip
    move = G4ThreeVector(0.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,+fFirstLayerZ/2.+fFirstLayerGapZ/2.) + fFirstLayerOffset;
    fAssemblyTISTAR->AddPlacedVolume(fFirstLayerLV,move,rotate);
    // Add first layer backward right strip
    move = G4ThreeVector(0.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,-fFirstLayerZ/2.-fFirstLayerGapZ/2.) + fFirstLayerOffset;
    fAssemblyTISTAR->AddPlacedVolume(fFirstLayerLV,move,rotate);

    // Build the second layer
    G4Box * secondLayerPV = new G4Box("secondLayerPV",fSecondLayerX/2.,fSecondLayerThickness/2.,fSecondLayerZ/2.);
    if(fSecondLayerLV == NULL) {
        fSecondLayerLV = new G4LogicalVolume(secondLayerPV,siliconMaterial,"secondLayerLV",0,0,0);
        fSecondLayerLV->SetVisAttributes(siliconVisAtt);
    }
    // Add second layer left strip
    move = G4ThreeVector(-fSecondLayerPCBUpperX/2.+fSecondLayerPCBLowerX/2.,+fSecondLayerThickness/2.+fSecondLayerDistFromBeam,0.) + fSecondLayerOffset;
    fAssemblyTISTAR->AddPlacedVolume(fSecondLayerLV,move,rotate);
    // Add second layer right strip
    move = G4ThreeVector(-fSecondLayerPCBUpperX/2.+fSecondLayerPCBLowerX/2.,-fSecondLayerThickness/2.-fSecondLayerDistFromBeam,0.) + fSecondLayerOffset;
    fAssemblyTISTAR->AddPlacedVolume(fSecondLayerLV,move,rotate);
    
    // Build the third layer
    G4Box * thirdLayerPV = new G4Box("thirdLayerPV",fThirdLayerX/2.,fThirdLayerThickness/2.,fThirdLayerZ/2.);
    if(fThirdLayerLV == NULL) {
        fThirdLayerLV = new G4LogicalVolume(thirdLayerPV,siliconMaterial,"thirdLayerLV",0,0,0);
        fThirdLayerLV->SetVisAttributes(siliconVisAtt);
    }
    // Add third left strip
    move = G4ThreeVector(-fThirdLayerPCBUpperX/2.+fThirdLayerPCBLowerX/2.,+fThirdLayerThickness/2.+fThirdLayerDistFromBeam,0.) + fThirdLayerOffset;
    fAssemblyTISTAR->AddPlacedVolume(fThirdLayerLV,move,rotate);
    // Add third right strip
    move = G4ThreeVector(-fThirdLayerPCBUpperX/2.+fThirdLayerPCBLowerX/2.,-fThirdLayerThickness/2.-fThirdLayerDistFromBeam,0.) + fThirdLayerOffset;
    fAssemblyTISTAR->AddPlacedVolume(fThirdLayerLV,move,rotate);

    return 1;
}

G4int DetectionSystemTISTAR::BuildPCBs()
{
    G4ThreeVector move, direction;
    G4RotationMatrix* rotate = new G4RotationMatrix();

    // Get the materials (should be built in the DetectorConstructionSuppressed)
    G4Material* PCBMaterial = G4Material::GetMaterial(fPCBMaterialName);
    if( !PCBMaterial ) {
        G4cout << " ----> Material " << fPCBMaterialName << " not found, cannot build! " << G4endl;
        return 0;
    }
    
    // Set up colours and other vis. attributes
    G4VisAttributes * PCBVisAtt = new G4VisAttributes(G4Colour::Green());
    PCBVisAtt->SetVisibility(true);

    // Build the first layer upper PCB
    if(fFirstLayerPCBUpperX > 0.) {
        G4Box * firstLayerPCBUpperXPV = new G4Box("firstLayerPCBUpperXPV",fFirstLayerPCBUpperX/2.,fPCBThickness/2.,fFirstLayerZ/2.);
        if(fFirstLayerPCBUpperXLV == NULL) {
            fFirstLayerPCBUpperXLV = new G4LogicalVolume(firstLayerPCBUpperXPV,PCBMaterial,"firstLayerPCBUpperXLV",0,0,0);
            fFirstLayerPCBUpperXLV->SetVisAttributes(PCBVisAtt);
        }
        // Add first layer upper forward left strip PCB
        move = G4ThreeVector(+fFirstLayerPCBUpperX/2.+fFirstLayerX/2.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,+fFirstLayerZ/2.+fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBUpperXLV,move,rotate);
        // Add first layer upper backward left strip PCB
        move = G4ThreeVector(+fFirstLayerPCBUpperX/2.+fFirstLayerX/2.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,-fFirstLayerZ/2.-fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBUpperXLV,move,rotate);
        // Add first layer upper forward right strip PCB
        move = G4ThreeVector(+fFirstLayerPCBUpperX/2.+fFirstLayerX/2.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,+fFirstLayerZ/2.+fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBUpperXLV,move,rotate);
        // Add first layer upper backward right strip PCB
        move = G4ThreeVector(+fFirstLayerPCBUpperX/2.+fFirstLayerX/2.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,-fFirstLayerZ/2.-fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBUpperXLV,move,rotate);
    }    
    // Build the first layer Lower PCB
    if(fFirstLayerPCBLowerX > 0.) {
        G4Box * firstLayerPCBLowerXPV = new G4Box("firstLayerPCBLowerXPV",fFirstLayerPCBLowerX/2.,fPCBThickness/2.,fFirstLayerZ/2.);
        if(fFirstLayerPCBLowerXLV == NULL) {
            fFirstLayerPCBLowerXLV = new G4LogicalVolume(firstLayerPCBLowerXPV,PCBMaterial,"firstLayerPCBLowerXLV",0,0,0);
            fFirstLayerPCBLowerXLV->SetVisAttributes(PCBVisAtt);
        }
        // Add first layer lower forward left strip PCB
        move = G4ThreeVector(-fFirstLayerPCBLowerX/2.-fFirstLayerX/2.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,+fFirstLayerZ/2.+fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBLowerXLV,move,rotate);
        // Add first layer lower backward left strip PCB
        move = G4ThreeVector(-fFirstLayerPCBLowerX/2.-fFirstLayerX/2.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,-fFirstLayerZ/2.-fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBLowerXLV,move,rotate);
        // Add first layer lower forward right strip PCB
        move = G4ThreeVector(-fFirstLayerPCBLowerX/2.-fFirstLayerX/2.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,+fFirstLayerZ/2.+fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBLowerXLV,move,rotate);
        // Add first layer lower backward right strip PCB
        move = G4ThreeVector(-fFirstLayerPCBLowerX/2.-fFirstLayerX/2.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,-fFirstLayerZ/2.-fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBLowerXLV,move,rotate);
    }
    
    // Build the first layer forward PCB
    if(fFirstLayerPCBForwardZ > 0.) {
        G4Box * firstLayerPCBForwardZPV = new G4Box("firstLayerPCBForwardZPV",fFirstLayerX/2.+fFirstLayerPCBUpperX/2.+fFirstLayerPCBLowerX/2.,fPCBThickness/2.,fFirstLayerPCBForwardZ/2.);
        if(fFirstLayerPCBForwardZLV == NULL) {
            fFirstLayerPCBForwardZLV = new G4LogicalVolume(firstLayerPCBForwardZPV,PCBMaterial,"firstLayerPCBForwardZLV",0,0,0);
            fFirstLayerPCBForwardZLV->SetVisAttributes(PCBVisAtt);
        }
        // Add first layer forward left strip PCB
        move = G4ThreeVector(+fFirstLayerPCBUpperX/2.-fFirstLayerPCBLowerX/2.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,+fFirstLayerPCBForwardZ/2.+fFirstLayerZ+fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBForwardZLV,move,rotate);
        // Add first layer forward right strip PCB
        move = G4ThreeVector(+fFirstLayerPCBUpperX/2.-fFirstLayerPCBLowerX/2.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,+fFirstLayerPCBForwardZ/2.+fFirstLayerZ+fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBForwardZLV,move,rotate);
    }
    // Build the first layer backward PCB
    if(fFirstLayerPCBBackwardZ > 0.) {
        G4Box * firstLayerPCBBackwardZPV = new G4Box("firstLayerPCBBackwardZPV",fFirstLayerX/2.+fFirstLayerPCBUpperX/2.+fFirstLayerPCBLowerX/2.,fPCBThickness/2.,fFirstLayerPCBBackwardZ/2.);
        if(fFirstLayerPCBBackwardZLV == NULL) {
            fFirstLayerPCBBackwardZLV = new G4LogicalVolume(firstLayerPCBBackwardZPV,PCBMaterial,"firstLayerPCBBackwardZLV",0,0,0);
            fFirstLayerPCBBackwardZLV->SetVisAttributes(PCBVisAtt);
        }
        // Add first layer forward left strip PCB
        move = G4ThreeVector(+fFirstLayerPCBUpperX/2.-fFirstLayerPCBLowerX/2.,+fFirstLayerThickness/2.+fFirstLayerDistFromBeam,-fFirstLayerPCBBackwardZ/2.-fFirstLayerZ-fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBBackwardZLV,move,rotate);
        // Add first layer forward right strip PCB
        move = G4ThreeVector(+fFirstLayerPCBUpperX/2.-fFirstLayerPCBLowerX/2.,-fFirstLayerThickness/2.-fFirstLayerDistFromBeam,-fFirstLayerPCBBackwardZ/2.-fFirstLayerZ-fFirstLayerGapZ/2.) + fFirstLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fFirstLayerPCBBackwardZLV,move,rotate);
    }

    // Build the Second layer upper PCB
    if(fSecondLayerPCBUpperX > 0.) {
        G4Box * secondLayerPCBUpperXPV = new G4Box("secondLayerPCBUpperXPV",fSecondLayerPCBUpperX/2.,fPCBThickness/2.,fSecondLayerZ/2.);
        if(fSecondLayerPCBUpperXLV == NULL) {
            fSecondLayerPCBUpperXLV = new G4LogicalVolume(secondLayerPCBUpperXPV,PCBMaterial,"secondLayerPCBUpperXLV",0,0,0);
            fSecondLayerPCBUpperXLV->SetVisAttributes(PCBVisAtt);
        }
        // Add Second layer upper forward left strip PCB
        move = G4ThreeVector(+fSecondLayerPCBLowerX/2.+fSecondLayerX/2.,+fSecondLayerThickness/2.+fSecondLayerDistFromBeam,0.) + fSecondLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fSecondLayerPCBUpperXLV,move,rotate);
        // Add Second layer upper backward right strip PCB
        move = G4ThreeVector(+fSecondLayerPCBLowerX/2.+fSecondLayerX/2.,-fSecondLayerThickness/2.-fSecondLayerDistFromBeam,0.) + fSecondLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fSecondLayerPCBUpperXLV,move,rotate);
    }    
    // Build the Second layer Lower PCB
    if(fSecondLayerPCBLowerX > 0.) {
        G4Box * secondLayerPCBLowerXPV = new G4Box("secondLayerPCBLowerXPV",fSecondLayerPCBLowerX/2.,fPCBThickness/2.,fSecondLayerZ/2.);
        if(fSecondLayerPCBLowerXLV == NULL) {
            fSecondLayerPCBLowerXLV = new G4LogicalVolume(secondLayerPCBLowerXPV,PCBMaterial,"secondLayerPCBLowerXLV",0,0,0);
            fSecondLayerPCBLowerXLV->SetVisAttributes(PCBVisAtt);
        }
        // Add Second layer lower forward left strip PCB
        move = G4ThreeVector(-fSecondLayerPCBUpperX/2.-fSecondLayerX/2.,+fSecondLayerThickness/2.+fSecondLayerDistFromBeam,0.) + fSecondLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fSecondLayerPCBLowerXLV,move,rotate);
        // Add Second layer lower backward right strip PCB
        move = G4ThreeVector(-fSecondLayerPCBUpperX/2.-fSecondLayerX/2.,-fSecondLayerThickness/2.-fSecondLayerDistFromBeam,0.) + fSecondLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fSecondLayerPCBLowerXLV,move,rotate);
    }

    // Build the second layer forward PCB
    if(fSecondLayerPCBForwardZ > 0.) {
        G4Box * secondLayerPCBForwardZPV = new G4Box("secondLayerPCBForwardZPV",fSecondLayerX/2.+fSecondLayerPCBUpperX/2.+fSecondLayerPCBLowerX/2.,fPCBThickness/2.,fSecondLayerPCBForwardZ/2.);
        if(fSecondLayerPCBForwardZLV == NULL) {
            fSecondLayerPCBForwardZLV = new G4LogicalVolume(secondLayerPCBForwardZPV,PCBMaterial,"secondLayerPCBForwardZLV",0,0,0);
            fSecondLayerPCBForwardZLV->SetVisAttributes(PCBVisAtt);
        }
        // Add Second layer forward left strip PCB
        move = G4ThreeVector(0.,+fSecondLayerThickness/2.+fSecondLayerDistFromBeam,+fSecondLayerPCBForwardZ/2.+fSecondLayerZ/2.) + fSecondLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fSecondLayerPCBForwardZLV,move,rotate);
        // Add Second layer forward right strip PCB
        move = G4ThreeVector(0.,-fSecondLayerThickness/2.-fSecondLayerDistFromBeam,+fSecondLayerPCBForwardZ/2.+fSecondLayerZ/2.) + fSecondLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fSecondLayerPCBForwardZLV,move,rotate);
    }
    // Build the second layer backward PCB
    if(fSecondLayerPCBBackwardZ > 0.) {
        G4Box * secondLayerPCBBackwardZPV = new G4Box("secondLayerPCBBackwardZPV",fSecondLayerX/2.+fSecondLayerPCBUpperX/2.+fSecondLayerPCBLowerX/2.,fPCBThickness/2.,fSecondLayerPCBBackwardZ/2.);
        if(fSecondLayerPCBBackwardZLV == NULL) {
            fSecondLayerPCBBackwardZLV = new G4LogicalVolume(secondLayerPCBBackwardZPV,PCBMaterial,"secondLayerPCBBackwardZLV",0,0,0);
            fSecondLayerPCBBackwardZLV->SetVisAttributes(PCBVisAtt);
        }
        // Add Second layer forward left strip PCB
        move = G4ThreeVector(0.,+fSecondLayerThickness/2.+fSecondLayerDistFromBeam,-fSecondLayerPCBBackwardZ/2.-fSecondLayerZ/2.) + fSecondLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fSecondLayerPCBBackwardZLV,move,rotate);
        // Add Second layer forward right strip PCB
        move = G4ThreeVector(0.,-fSecondLayerThickness/2.-fSecondLayerDistFromBeam,-fSecondLayerPCBBackwardZ/2.-fSecondLayerZ/2.) + fSecondLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fSecondLayerPCBBackwardZLV,move,rotate);
    }
    
    // Build the Third layer upper PCB
    if(fThirdLayerPCBUpperX > 0.) {
        G4Box * thirdLayerPCBUpperXPV = new G4Box("thirdLayerPCBUpperXPV",fThirdLayerPCBUpperX/2.,fPCBThickness/2.,fThirdLayerZ/2.);
        if(fThirdLayerPCBUpperXLV == NULL) {
            fThirdLayerPCBUpperXLV = new G4LogicalVolume(thirdLayerPCBUpperXPV,PCBMaterial,"thirdLayerPCBUpperXLV",0,0,0);
            fThirdLayerPCBUpperXLV->SetVisAttributes(PCBVisAtt);
        }
        // Add Third layer upper forward left strip PCB
        move = G4ThreeVector(+fThirdLayerPCBLowerX/2.+fThirdLayerX/2.,+fThirdLayerThickness/2.+fThirdLayerDistFromBeam,0.) + fThirdLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fThirdLayerPCBUpperXLV,move,rotate);
        // Add Third layer upper backward right strip PCB
        move = G4ThreeVector(+fThirdLayerPCBLowerX/2.+fThirdLayerX/2.,-fThirdLayerThickness/2.-fThirdLayerDistFromBeam,0.) + fThirdLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fThirdLayerPCBUpperXLV,move,rotate);
    }    
    // Build the Third layer Lower PCB
    if(fThirdLayerPCBLowerX > 0.) {
        G4Box * thirdLayerPCBLowerXPV = new G4Box("thirdLayerPCBLowerXPV",fThirdLayerPCBLowerX/2.,fPCBThickness/2.,fThirdLayerZ/2.);
        if(fThirdLayerPCBLowerXLV == NULL) {
            fThirdLayerPCBLowerXLV = new G4LogicalVolume(thirdLayerPCBLowerXPV,PCBMaterial,"thirdLayerPCBLowerXLV",0,0,0);
            fThirdLayerPCBLowerXLV->SetVisAttributes(PCBVisAtt);
        }
        // Add Third layer lower forward left strip PCB
        move = G4ThreeVector(-fThirdLayerPCBUpperX/2.-fThirdLayerX/2.,+fThirdLayerThickness/2.+fThirdLayerDistFromBeam,0.) + fThirdLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fThirdLayerPCBLowerXLV,move,rotate);
        // Add Third layer lower backward right strip PCB
        move = G4ThreeVector(-fThirdLayerPCBUpperX/2.-fThirdLayerX/2.,-fThirdLayerThickness/2.-fThirdLayerDistFromBeam,0.) + fThirdLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fThirdLayerPCBLowerXLV,move,rotate);
    }
    
    // Build the third layer forward PCB
    if(fThirdLayerPCBForwardZ > 0.) {
        G4Box * thirdLayerPCBForwardZPV = new G4Box("thirdLayerPCBForwardZPV",fThirdLayerX/2.+fThirdLayerPCBUpperX/2.+fThirdLayerPCBLowerX/2.,fPCBThickness/2.,fThirdLayerPCBForwardZ/2.);
        if(fThirdLayerPCBForwardZLV == NULL) {
            fThirdLayerPCBForwardZLV = new G4LogicalVolume(thirdLayerPCBForwardZPV,PCBMaterial,"thirdLayerPCBForwardZLV",0,0,0);
            fThirdLayerPCBForwardZLV->SetVisAttributes(PCBVisAtt);
        }
        // Add Third layer forward left strip PCB
        move = G4ThreeVector(0.,+fThirdLayerThickness/2.+fThirdLayerDistFromBeam,+fThirdLayerPCBForwardZ/2.+fThirdLayerZ/2.) + fThirdLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fThirdLayerPCBForwardZLV,move,rotate);
        // Add Third layer forward right strip PCB
        move = G4ThreeVector(0.,-fThirdLayerThickness/2.-fThirdLayerDistFromBeam,+fThirdLayerPCBForwardZ/2.+fThirdLayerZ/2.) + fThirdLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fThirdLayerPCBForwardZLV,move,rotate);
    }
    // Build the third layer backward PCB
    if(fThirdLayerPCBBackwardZ > 0.) {
        G4Box * thirdLayerPCBBackwardZPV = new G4Box("thirdLayerPCBBackwardZPV",fThirdLayerX/2.+fThirdLayerPCBUpperX/2.+fThirdLayerPCBLowerX/2.,fPCBThickness/2.,fThirdLayerPCBBackwardZ/2.);
        if(fThirdLayerPCBBackwardZLV == NULL) {
            fThirdLayerPCBBackwardZLV = new G4LogicalVolume(thirdLayerPCBBackwardZPV,PCBMaterial,"thirdLayerPCBBackwardZLV",0,0,0);
            fThirdLayerPCBBackwardZLV->SetVisAttributes(PCBVisAtt);
        }
        // Add Third layer forward left strip PCB
        move = G4ThreeVector(0.,+fThirdLayerThickness/2.+fThirdLayerDistFromBeam,-fThirdLayerPCBBackwardZ/2.-fThirdLayerZ/2.) + fThirdLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fThirdLayerPCBBackwardZLV,move,rotate);
        // Add Third layer forward right strip PCB
        move = G4ThreeVector(0.,-fThirdLayerThickness/2.-fThirdLayerDistFromBeam,-fThirdLayerPCBBackwardZ/2.-fThirdLayerZ/2.) + fThirdLayerOffset;
        fAssemblyTISTAR->AddPlacedVolume(fThirdLayerPCBBackwardZLV,move,rotate);
    }

    return 1;
}

G4int DetectionSystemTISTAR::BuildLayer(G4int layer) 
{
    if(layer<0||layer>=fNLayers) { 
        G4cout << " ---> This layer doesn't exist! " << G4endl;
        return 0;
    }
    if(!fDimensionsSet.at(layer).at(0) || !fDimensionsSet.at(layer).at(1)) {
        G4cout << " ---> Dimensions not set!" << G4endl;
        return 0;
    }

    G4String name;
    G4RotationMatrix * rotate = new G4RotationMatrix();
    G4ThreeVector move;    

    // Get materials
    G4Material* si_material = G4Material::GetMaterial(fSiliconMaterialName);
    if( !si_material ) { G4cout << " ---> Material " << fSiliconMaterialName << " not found, cannot build! " << G4endl; return 0; }
    G4Material* pcb_material = G4Material::GetMaterial(fPCBMaterialName);
    if( !pcb_material ) { G4cout << " ---> Material " << fPCBMaterialName << " not found, cannot build! " << G4endl; return 0; }
   
    // Set up colours and other vis. attributes
    G4VisAttributes * si_vis_att = new G4VisAttributes(G4Colour::Cyan());
    si_vis_att->SetVisibility(true);
    G4VisAttributes * pcb_vis_att = new G4VisAttributes(G4Colour::Green());
    pcb_vis_att->SetVisibility(true);

    // Get/calculate dimensions
    G4ThreeVector si_dim = fDimensionsSiLayers.at(layer);
    G4ThreeVector pcb_dim = fDimensionsPCBLayers.at(layer);
    G4ThreeVector offset = fOffsetLayers.at(layer); // offset of Si layer from top left corner of PCB layer
    
    G4ThreeVector cut_extra = G4ThreeVector(0.,1.*m,0.);
    if(offset.x()<=0.) cut_extra += G4ThreeVector(1.*m,0.,0.);
    if(offset.z()<=0.) cut_extra += G4ThreeVector(0.,0.,1.*m);
    
    move = G4ThreeVector( -si_dim.x()/2. +pcb_dim.x()/2. -offset.x() +cut_extra.x(),  // x
                          0.,                                                         // y  
                          +si_dim.z()/2. -pcb_dim.z()/2. +offset.z() -cut_extra.z()); // z

    // Physical volumes
    
    name  = "TISTARSiLayer" + std::to_string(layer) + "PV";
    G4Box * si_PV = new G4Box(name,si_dim.x()/2.,si_dim.y()/2.,si_dim.z()/2.);
    
    name = "TISTARPCBLayerPreCut" + std::to_string(layer) + "PV";
    G4Box * pcb_precut_PV = new G4Box(name,pcb_dim.x()/2.,pcb_dim.y()/2.,pcb_dim.z()/2.);    
    
    name = "TISTARPCBLayer" + std::to_string(layer) + "CuttingPV";
    G4Box * pcb_cutting_PV = new G4Box(name,si_dim.x()/2.+cut_extra.x(),pcb_dim.y()/2.+cut_extra.y(),si_dim.z()/2.+cut_extra.z());
    
    name = "TISTARPCBLayer" + std::to_string(layer) + "PV";
    G4SubtractionSolid * pcb_PV = new G4SubtractionSolid(name,pcb_precut_PV,pcb_cutting_PV,rotate,move);
    
    // Logical volumes
    if(fLogicalSiLayers.at(layer) == NULL) {
        name = "TISTARSiLayer" + std::to_string(layer) + "LV";
        fLogicalSiLayers.at(layer) = new G4LogicalVolume(si_PV,si_material,name,0,0,0);
        fLogicalSiLayers.at(layer)->SetVisAttributes(si_vis_att);
    }
    else { G4cout << " ---> Already built this Si layer - " << layer << G4endl; return 0; }
    
    if(fLogicalPCBLayers.at(layer) == NULL) {
        name = "TISTARPCBLayer" + std::to_string(layer) + "LV";
        fLogicalPCBLayers.at(layer) = new G4LogicalVolume(pcb_PV,pcb_material,name,0,0,0);
        fLogicalPCBLayers.at(layer)->SetVisAttributes(pcb_vis_att);
    }
    else { G4cout << " ---> Already built this PCB layer - " << layer << G4endl; return 0; }
    
    // Add to assembly volume
    if(fAssemblyLayers.at(layer) == NULL) {
        fAssemblyLayers.at(layer) = new G4AssemblyVolume();
        // Si layer
        move = G4ThreeVector(0.,0.,0.);
        fAssemblyLayers.at(layer)->AddPlacedVolume(fLogicalSiLayers.at(layer),move,rotate);
        // PCB layer
        move = G4ThreeVector( -pcb_dim.x()/2. + si_dim.x()/2. + offset.x(),  // x
                              0.,                                            // y
                              +pcb_dim.z()/2. - si_dim.z()/2. - offset.z()); // z
        fAssemblyLayers.at(layer)->AddPlacedVolume(fLogicalPCBLayers.at(layer),move,rotate);
    }

    return 1;
}

