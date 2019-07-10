#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
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
#include "TRexSettings.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <string>

DetectionSystemTISTAR::DetectionSystemTISTAR() :
    fAssemblyLayer(NULL),
    fLogicalSiLayer(NULL),
    fLogicalPCBLayer(NULL),
    fLogicalGasTarget(NULL),
    fLogicalGasTargetMylar(NULL),
    fLogicalGasTargetBeWindow(NULL),
    fSiDimensionsSet(false),
    fPCBDimensionsSet(false)
{
    fSiliconMaterialName = "Silicon";
    fPCBMaterialName = "Epoxy-Resin";

    fSiDimensions = G4ThreeVector(0.,0.,0.);
    fPCBDimensions = G4ThreeVector(0.,0.,0.);
    fSiOffsetInPCB = G4ThreeVector(0.,0.,0.);
    fPositionOffset = G4ThreeVector(0.,0.,0.);

    fGasTargetRadius = TRexSettings::Get()->GetTargetDiameter()/2.;
    fGasTargetLength = TRexSettings::Get()->GetGasTargetLength();
    fGasTargetDensity = TRexSettings::Get()->GetTargetMaterialDensity();
    fGasTargetPressure = TRexSettings::Get()->GetTargetPressure();
    fGasTargetMaterialName = TRexSettings::Get()->GetTargetMaterialName();

    fGasTargetMylarThickness = TRexSettings::Get()->GetTargetMylarThickness();
    fGasTargetMylarMaterialName = "Mylar";

    fGasTargetBeWindowThickness = TRexSettings::Get()->GetTargetBeWindowThickness();
    fGasTargetBeWindowMaterialName = "Beryllium";
    
    fDetectorNumber = 10;

    fVacuumChamberShape = "cylinder";
    fVacuumChamberMaterialName = "Vacuum";
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

    fAssemblyLayer->MakeImprint(expHallLog, move, rotate);

    return 1;
}

G4int DetectionSystemTISTAR::PlaceDetector(G4ThreeVector move, G4ThreeVector rotate, G4LogicalVolume* expHallLog)
{
    G4RotationMatrix * rotation = new G4RotationMatrix;
    rotation->rotateX(rotate.x()*M_PI/180.);
    rotation->rotateY(rotate.y()*M_PI/180.);
    rotation->rotateZ(rotate.z()*M_PI/180.);
    fAssemblyLayer->MakeImprint(expHallLog, move, rotation);

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
    G4ThreeVector offset = fSiOffsetInPCB; // offset of Si layer from top left corner of PCB layer

    G4ThreeVector cut_extra = G4ThreeVector(0.,1.*m,0.);
    if(offset.x()<=0.) cut_extra += G4ThreeVector(1.*m,0.,0.);
    if(offset.z()<=0.) cut_extra += G4ThreeVector(0.,0.,1.*m);

    move = G4ThreeVector( -si_dim.x()/2. +pcb_dim.x()/2. -offset.x() +cut_extra.x(),  // x
                          0.,                                                         // y  
                          +si_dim.z()/2. -pcb_dim.z()/2. +offset.z() -cut_extra.z()); // z

    // Physical volumes

    name  = "TISTARSiLayerPV";
    G4Box * si_PV = new G4Box(name,si_dim.x()/2.,si_dim.y()/2.,si_dim.z()/2.);

    G4Box * pcb_precut_PV = NULL;
    G4Box * pcb_cutting_PV = NULL;
    G4SubtractionSolid * pcb_PV = NULL;
    if(fPCBDimensionsSet) {
        name = "TISTARPCBLayerPreCutPV";
        pcb_precut_PV = new G4Box(name,pcb_dim.x()/2.,pcb_dim.y()/2.,pcb_dim.z()/2.);    
        name = "TISTARPCBLayerCuttingPV";
        pcb_cutting_PV = new G4Box(name,si_dim.x()/2.+cut_extra.x(),pcb_dim.y()/2.+cut_extra.y(),si_dim.z()/2.+cut_extra.z());
        name = "TISTARPCBLayerPV";
        pcb_PV = new G4SubtractionSolid(name,pcb_precut_PV,pcb_cutting_PV,rotate,move);
    }

    // Logical volumes
    if(fLogicalSiLayer == NULL) {
        name = "TISTARSiLayerLV";
        fLogicalSiLayer = new G4LogicalVolume(si_PV,si_material,name,0,0,0);
        fLogicalSiLayer->SetVisAttributes(si_vis_att);
    }

    if(fPCBDimensionsSet) {
        if(fLogicalPCBLayer == NULL) {
            name = "TISTARPCBLayerLV";
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

    G4String name;

    // first strip
    move = G4ThreeVector(   +dist_from_beam,        // x
                            0.,                     // y
                            0.);                    // z
    if(fPCBDimensionsSet && !si_centered) move += G4ThreeVector(0., -fSiOffsetInPCB.x() -fSiDimensions.x()/2. +fPCBDimensions.x()/2., 0.);
    rotate = G4ThreeVector( 0.,                     // x-rotation
                            0.,                     // y-rotation
                            +90.);                  // z-rotation
    name = "TISTARSiLayerLV_" + std::to_string(fDetectorNumber);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);

    // second strip
    move = G4ThreeVector(   -dist_from_beam,        // x
                            0.,                     // y
                            0.);                    // z
    if(fPCBDimensionsSet && !si_centered) move += G4ThreeVector(0., +fSiOffsetInPCB.x() +fSiDimensions.x()/2. -fPCBDimensions.x()/2., 0.);
    rotate = G4ThreeVector( 0,                      // x-rotation
                            +180.,                  // y-rotation
                            +90.);                  // z-rotation
    name = "TISTARSiLayerLV_" + std::to_string(fDetectorNumber+1);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);

    return 1;
}

G4int DetectionSystemTISTAR::Add4StripLayer(G4double dist_from_beam, G4double gap_z, G4LogicalVolume* expHallLog) 
{
    G4ThreeVector move;
    G4ThreeVector rotate;

    G4String name;

    //fPositionOffset = G4ThreeVector(0., 10.0*mm*(-5.0*mm -100.0/2.*mm +140.0/2.*mm)/(30.0*mm), 0.);

    // first strip
    move = G4ThreeVector(   +dist_from_beam,                // x
                            0.,                             // y
                            +fSiDimensions.z()/2. +gap_z);  // z
    move += fPositionOffset;
    rotate = G4ThreeVector( 0.,                     // x-rotation
                            0.,                     // y-rotation
                            +90.);                  // z-rotation
    name = "TISTARSiLayerLV_" + std::to_string(fDetectorNumber);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);

    // second strip
    move = G4ThreeVector(   +dist_from_beam,                // x
                            0.,                             // y
                            -fSiDimensions.z()/2. -gap_z);  // z
    move += fPositionOffset;
    rotate = G4ThreeVector( 0,                      // x-rotation
                            +180.,                  // y-rotation
                            +90.);                  // z-rotation
    name = "TISTARSiLayerLV_" + std::to_string(fDetectorNumber+1);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);

    // third strip
    move = G4ThreeVector(   -dist_from_beam,                // x
                            0.,                             // y
                            +fSiDimensions.z()/2. +gap_z);  // z
    move -= fPositionOffset;
    rotate = G4ThreeVector( 0.,                     // x-rotation
                            0.,                     // y-rotation
                            -90.);                  // z-rotation
    name = "TISTARSiLayerLV_" + std::to_string(fDetectorNumber+2);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);

    // fourth strip
    move = G4ThreeVector(   -dist_from_beam,                // x
                            0.,                             // y
                            -fSiDimensions.z()/2. -gap_z);  // z
    move -= fPositionOffset;
    rotate = G4ThreeVector( 0,                      // x-rotation
                            +180.,                  // y-rotation
                            -90.);                  // z-rotation
    name = "TISTARSiLayerLV_" + std::to_string(fDetectorNumber+3);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);


    return 1;
}

G4int DetectionSystemTISTAR::AddGasTarget(G4LogicalVolume* expHallLog)
{
    G4ThreeVector move;
    G4RotationMatrix * rotate = NULL;

    // Make the target materials
    G4Material * deuterium_gas_material = new G4Material("Deuterium_Gas", 1, 2.014*g/mole, fGasTargetDensity, kStateGas, 293.15*kelvin, fGasTargetPressure);
    G4Material * mylar_material = G4Material::GetMaterial(fGasTargetMylarMaterialName);    
    G4Material * be_material = G4Material::GetMaterial(fGasTargetBeWindowMaterialName);    
    std::cout<<"TI-STAR gas target built w/ density = "<<G4BestUnit(fGasTargetDensity,"Volumic Mass")<<" and pressure = "<<G4BestUnit(fGasTargetPressure,"Pressure")<<std::endl;

    // Set up colours and other vis. attributes
    G4VisAttributes * gas_vis_att = new G4VisAttributes(G4Colour::Yellow());
    gas_vis_att->SetVisibility(true);

    G4VisAttributes * mylar_vis_att = new G4VisAttributes(G4Colour::Magenta());
    mylar_vis_att->SetVisibility(true);

    G4VisAttributes * be_vis_att = new G4VisAttributes(G4Colour::Blue());
    be_vis_att->SetVisibility(true);
 
    // Gas    
    // Build the object volume
    G4Tubs * gas_target = new G4Tubs("Gas_target", 0.*cm, fGasTargetRadius, fGasTargetLength/2., 0., 2.*M_PI);
    // Logical volumes
    if(fLogicalGasTarget == NULL) {
        fLogicalGasTarget = new G4LogicalVolume(gas_target,deuterium_gas_material,"Gas_target_LV",0,0,0);
        fLogicalGasTarget->SetVisAttributes(gas_vis_att);
    }
    // Placement
    move = G4ThreeVector(0.,0.,0.);
    rotate = new G4RotationMatrix;
    G4VPhysicalVolume * gas_target_PV = new G4PVPlacement(rotate, move, fLogicalGasTarget, "Gas_target_PV", expHallLog, 0, 0, 0);

    // Mylar foil
    // Build the object volume
    G4Tubs * mylar = new G4Tubs("Gas_target_mylar", fGasTargetRadius, fGasTargetRadius+fGasTargetMylarThickness, fGasTargetLength/2., 0., 2.*M_PI);
    // Logical volumes
    if(fLogicalGasTargetMylar == NULL) {
        fLogicalGasTargetMylar = new G4LogicalVolume(mylar,mylar_material,"Gas_target_mylar_LV",0,0,0);
        fLogicalGasTargetMylar->SetVisAttributes(mylar_vis_att);
    }
    // Placement
    move = G4ThreeVector(0.,0.,0.);
    rotate = new G4RotationMatrix;
    G4VPhysicalVolume * gas_target_mylar_PV = new G4PVPlacement(rotate, move, fLogicalGasTargetMylar, "Gas_target_mylar_PV", expHallLog, 0, 0, 0);

    // Be window
    G4Tubs * be_window = new G4Tubs("Gas_target_Be_window", 0.*cm, fGasTargetRadius+fGasTargetMylarThickness, fGasTargetBeWindowThickness/2., 0., 2.*M_PI);
    // Logical volumes
    if(fLogicalGasTargetBeWindow == NULL) {
        fLogicalGasTargetBeWindow = new G4LogicalVolume(be_window, be_material,"Gas_target_Be_Window_LV",0,0,0);
        fLogicalGasTargetBeWindow->SetVisAttributes(be_vis_att);
    }
    // Placement
    // front
    move = G4ThreeVector(0., 0., +fGasTargetBeWindowThickness/2.+fGasTargetLength/2.);
    rotate = new G4RotationMatrix;
    G4VPhysicalVolume * gas_target_mylar_PV1 = new G4PVPlacement(rotate, move, fLogicalGasTargetBeWindow, "Gas_target_mylar_PV1", expHallLog, 0, 0, 0);
    // back
    move = G4ThreeVector(0., 0., -fGasTargetBeWindowThickness/2.-fGasTargetLength/2.);
    rotate = new G4RotationMatrix;
    G4VPhysicalVolume * gas_target_mylar_PV2 = new G4PVPlacement(rotate, move, fLogicalGasTargetBeWindow, "Gas_target_mylar_PV2", expHallLog, 0, 0, 0);

    return 1;
}

void DetectionSystemTISTAR::SetVacuumChamberShape(G4String shape) 
{
    if(shape=="box" || shape=="cylinder") {
        fVacuumChamberShape = shape;
    } else {
        G4cout << " ---> TI-STAR vacuum chamber shape \"" << shape << "\" unknown!" << G4endl;
    }
}

G4int DetectionSystemTISTAR::AddVacuumChamber(G4LogicalVolume* expHallLog, G4LogicalVolume *& vacuumChamberLog) 
{
    G4ThreeVector move;
    G4RotationMatrix * rotate = NULL;

    // Make the target materials
    G4Material * vacuum_material = G4Material::GetMaterial(fVacuumChamberMaterialName);

    // Set up colours and other vis. attributes
    G4VisAttributes * vacuum_vis_att = new G4VisAttributes(G4Colour::Red());
    vacuum_vis_att->SetForceWireframe(true);
    vacuum_vis_att->SetVisibility(true);

    // Build the object solid volume
    G4VSolid * vacuum_chamber_SV = NULL;
    if(fVacuumChamberShape == "box") {
        vacuum_chamber_SV = new G4Box("vacuum_chamber_solid", fVacuumChamberBoxDimensions.x(), fVacuumChamberBoxDimensions.y(), fVacuumChamberBoxDimensions.z());
    } else if(fVacuumChamberShape == "cylinder") {
        vacuum_chamber_SV = new G4Tubs("vacuum_chamber_solid", 0., fVacuumChamberCylinderRadius, fVacuumChamberCylinderZ/2., 0, 2.*M_PI);
    } else {
        G4cout << " ---> Unknown vacuum chamber shape \"" << fVacuumChamberShape << "\", cannot build!" << G4endl;
        return 2;
    }
    
    // Logical volumes
    vacuumChamberLog = new G4LogicalVolume(vacuum_chamber_SV, vacuum_material, "vacuum_chamber_LV", 0, 0, 0);
    vacuumChamberLog->SetVisAttributes(vacuum_vis_att);

    // Placement
    move = G4ThreeVector(0.,0.,0.);
    rotate = new G4RotationMatrix;
    G4VPhysicalVolume * vacuum_chamber_PV = new G4PVPlacement(rotate, move, vacuumChamberLog, "vacuum_chamber_PV", expHallLog, 0, 0, 0);    

    return 1;
}

