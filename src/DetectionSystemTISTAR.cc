#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Box.hh"
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

DetectionSystemTISTAR::DetectionSystemTISTAR(G4int layer_number) :
    fLayerNumber(layer_number),
    fAssemblyLayer(NULL),
    fLogicalSiLayer(NULL),
    fLogicalPCBLayer(NULL),
    fSiDimensionsSet(false),
    fPCBDimensionsSet(false)
{
    fSiliconMaterialName = "Silicon";
    fPCBMaterialName = "Epoxy-Resin";

    fSiDimensions = G4ThreeVector(0.,0.,0.);
    fPCBDimensions = G4ThreeVector(0.,0.,0.);
    fOffset = G4ThreeVector(0.,0.,0.);

}

DetectionSystemTISTAR::~DetectionSystemTISTAR() 
{
    delete fLogicalSiLayer;
    delete fLogicalPCBLayer;
    delete fAssemblyLayer;
}

G4int DetectionSystemTISTAR::Build() 
{
    BuildLayer();

    return 1;
}

G4int DetectionSystemTISTAR::PlaceDetector(G4LogicalVolume* expHallLog) 
{
    G4RotationMatrix * rotate = new G4RotationMatrix;
    G4ThreeVector move = G4ThreeVector(0., 0., 0.);

    fAssemblyLayer->MakeImprint(expHallLog, move, rotate, fLayerNumber);

    return 1;
}

G4int DetectionSystemTISTAR::PlaceDetector(G4ThreeVector move, G4ThreeVector rotate, G4LogicalVolume* expHallLog)
{
    G4RotationMatrix * rotation = new G4RotationMatrix;
    rotation->rotateX(rotate.x()*M_PI/180.);
    rotation->rotateY(rotate.y()*M_PI/180.);
    rotation->rotateZ(rotate.z()*M_PI/180.);
    fAssemblyLayer->MakeImprint(expHallLog, move, rotation, fLayerNumber);

    return 1;
}

G4int DetectionSystemTISTAR::BuildLayer() 
{

    if(!fSiDimensionsSet) { G4cout << " ---> Silicon dimensions not set!" << G4endl; return 0; }

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
    G4ThreeVector si_dim = fSiDimensions;
    G4ThreeVector pcb_dim = fPCBDimensions;
    G4ThreeVector offset = fOffset; // offset of Si layer from top left corner of PCB layer

    G4ThreeVector cut_extra = G4ThreeVector(0.,1.*m,0.);
    if(offset.x()<=0.) cut_extra += G4ThreeVector(1.*m,0.,0.);
    if(offset.z()<=0.) cut_extra += G4ThreeVector(0.,0.,1.*m);

    move = G4ThreeVector( -si_dim.x()/2. +pcb_dim.x()/2. -offset.x() +cut_extra.x(),  // x
                          0.,                                                         // y  
                          +si_dim.z()/2. -pcb_dim.z()/2. +offset.z() -cut_extra.z()); // z

    // Physical volumes

    name  = "TISTARSiLayer" + std::to_string(fLayerNumber) + "PV";
    G4Box * si_PV = new G4Box(name,si_dim.x()/2.,si_dim.y()/2.,si_dim.z()/2.);

    G4Box * pcb_precut_PV = NULL;
    G4Box * pcb_cutting_PV = NULL;
    G4SubtractionSolid * pcb_PV = NULL;
    if(fPCBDimensionsSet) {
        name = "TISTARPCBLayerPreCut" + std::to_string(fLayerNumber) + "PV";
        pcb_precut_PV = new G4Box(name,pcb_dim.x()/2.,pcb_dim.y()/2.,pcb_dim.z()/2.);    
        name = "TISTARPCBLayer" + std::to_string(fLayerNumber) + "CuttingPV";
        pcb_cutting_PV = new G4Box(name,si_dim.x()/2.+cut_extra.x(),pcb_dim.y()/2.+cut_extra.y(),si_dim.z()/2.+cut_extra.z());
        name = "TISTARPCBLayer" + std::to_string(fLayerNumber) + "PV";
        pcb_PV = new G4SubtractionSolid(name,pcb_precut_PV,pcb_cutting_PV,rotate,move);
    }

    // Logical volumes
    if(fLogicalSiLayer == NULL) {
        name = "TISTARSiLayer" + std::to_string(fLayerNumber) + "LV";
        fLogicalSiLayer = new G4LogicalVolume(si_PV,si_material,name,0,0,0);
        fLogicalSiLayer->SetVisAttributes(si_vis_att);
    }

    if(fPCBDimensionsSet) {
        if(fLogicalPCBLayer == NULL) {
            name = "TISTARPCBLayer" + std::to_string(fLayerNumber) + "LV";
            fLogicalPCBLayer = new G4LogicalVolume(pcb_PV,pcb_material,name,0,0,0);
            fLogicalPCBLayer->SetVisAttributes(pcb_vis_att);
        }
    }    

    // Add to assembly volume
    if(fAssemblyLayer == NULL) {
        fAssemblyLayer = new G4AssemblyVolume();

        // Si layer
        move = G4ThreeVector(0.,0.,0.);
        fAssemblyLayer->AddPlacedVolume(fLogicalSiLayer,move,rotate);

        // PCB layer
        if(fPCBDimensionsSet) {
            move = G4ThreeVector( -pcb_dim.x()/2. + si_dim.x()/2. + offset.x(),  // x
                                  0.,                                            // y
                                  +pcb_dim.z()/2. - si_dim.z()/2. - offset.z()); // z
            fAssemblyLayer->AddPlacedVolume(fLogicalPCBLayer,move,rotate);
        }
    }

    return 1;
}

G4int DetectionSystemTISTAR::Add2StripLayer(G4double dist_from_beam, G4bool si_centered, G4LogicalVolume* expHallLog) 
{
    G4ThreeVector move;
    G4ThreeVector rotate;

    // first strip
    move = G4ThreeVector(   +dist_from_beam,        // x
                            0.,                     // y
                            0.);                    // z
    if(fPCBDimensionsSet && !si_centered) move += G4ThreeVector(0., -fOffset.x() -fSiDimensions.x()/2. +fPCBDimensions.x()/2., 0.);
    rotate = G4ThreeVector( 0.,                     // x-rotation
                            0.,                     // y-rotation
                            +90.);                  // z-rotation
    PlaceDetector(move, rotate, expHallLog);
    fLayerNumber++;

    // second strip
    move = G4ThreeVector(   -dist_from_beam,        // x
                            0.,                     // y
                            0.);                    // z
    if(fPCBDimensionsSet && !si_centered) move += G4ThreeVector(0., +fOffset.x() +fSiDimensions.x()/2. -fPCBDimensions.x()/2., 0.);
    rotate = G4ThreeVector( 0,                      // x-rotation
                            +180.,                  // y-rotation
                            +90.);                  // z-rotation
    PlaceDetector(move, rotate, expHallLog);

    return 1;
}

G4int DetectionSystemTISTAR::Add4StripLayer(G4double dist_from_beam, G4double gap_z, G4LogicalVolume* expHallLog) 
{
    G4ThreeVector move;
    G4ThreeVector rotate;

    // first strip
    move = G4ThreeVector(   +dist_from_beam,                // x
                            0.,                             // y
                            +fSiDimensions.z()/2. +gap_z);  // z
    rotate = G4ThreeVector( 0.,                     // x-rotation
                            0.,                     // y-rotation
                            +90.);                  // z-rotation
    PlaceDetector(move, rotate, expHallLog);
    fLayerNumber++;

    // second strip
    move = G4ThreeVector(   +dist_from_beam,                // x
                            0.,                             // y
                            -fSiDimensions.z()/2. -gap_z);  // z
    rotate = G4ThreeVector( 0,                      // x-rotation
                            +180.,                  // y-rotation
                            +90.);                  // z-rotation
    PlaceDetector(move, rotate, expHallLog);
    fLayerNumber++;

    // third strip
    move = G4ThreeVector(   -dist_from_beam,                // x
                            0.,                             // y
                            +fSiDimensions.z()/2. +gap_z);  // z
    rotate = G4ThreeVector( 0.,                     // x-rotation
                            0.,                     // y-rotation
                            -90.);                  // z-rotation
    PlaceDetector(move, rotate, expHallLog);
    fLayerNumber++;

    // fourth strip
    move = G4ThreeVector(   -dist_from_beam,                // x
                            0.,                             // y
                            -fSiDimensions.z()/2. -gap_z);  // z
    rotate = G4ThreeVector( 0,                      // x-rotation
                            +180.,                  // y-rotation
                            -90.);                  // z-rotation
    PlaceDetector(move, rotate, expHallLog);


    return 1;
}
