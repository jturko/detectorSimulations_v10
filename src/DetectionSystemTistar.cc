#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
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

#include "DetectionSystemTistar.hh"
#include "TistarSettings.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <string>

#include "G4NistManager.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"

DetectionSystemTistar::DetectionSystemTistar() :
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

    fGasTargetRadius = TistarSettings::Get()->GetTargetDiameter()/2.;
    fGasTargetLength = TistarSettings::Get()->GetGasTargetLength();
    fGasTargetDensity = TistarSettings::Get()->GetTargetMaterialDensity();
    fGasTargetPressure = TistarSettings::Get()->GetTargetPressure();
    fGasTargetMaterialName = TistarSettings::Get()->GetTargetMaterialName();

    fGasTargetMylarThickness = TistarSettings::Get()->GetTargetMylarThickness();
    fGasTargetMylarMaterialName = "Mylar";

    fGasTargetBeWindowThickness = TistarSettings::Get()->GetTargetBeWindowThickness();
    fGasTargetBeWindowMaterialName = "Beryllium";
    
    fDetectorNumber = 10;

    fVacuumChamberShape = "cylinder";
    fVacuumChamberMaterialName = "Vacuum";

    fVacuumChamberExteriorMaterialName = "G4_Al";
    fVacuumChamberExteriorThickness = 1.0*mm;
}

DetectionSystemTistar::~DetectionSystemTistar() 
{
    delete fLogicalSiLayer;
    delete fLogicalPCBLayer;
    delete fAssemblyLayer;
}

G4int DetectionSystemTistar::Build() 
{
    BuildLayer();

    return 1;
}

G4int DetectionSystemTistar::PlaceDetector(G4LogicalVolume* expHallLog) 
{
    G4RotationMatrix * rotate = new G4RotationMatrix;
    G4ThreeVector move = G4ThreeVector(0., 0., 0.);

    fAssemblyLayer->MakeImprint(expHallLog, move, rotate);

    return 1;
}

G4int DetectionSystemTistar::PlaceDetector(G4ThreeVector move, G4ThreeVector rotate, G4LogicalVolume* expHallLog)
{
    G4RotationMatrix * rotation = new G4RotationMatrix;
    rotation->rotateX(rotate.x()*M_PI/180.);
    rotation->rotateY(rotate.y()*M_PI/180.);
    rotation->rotateZ(rotate.z()*M_PI/180.);
    fAssemblyLayer->MakeImprint(expHallLog, move, rotation);

    return 1;
}

G4int DetectionSystemTistar::BuildLayer() 
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

    G4ThreeVector cut_extra = G4ThreeVector(1.*m,0.,0.);
    if(offset.y()<=0.) cut_extra += G4ThreeVector(0.,1.*m,0.);
    if(offset.z()<=0.) cut_extra += G4ThreeVector(0.,0.,1.*m);

    move = G4ThreeVector( 0.,                                                         // x
                          -si_dim.y()/2. +pcb_dim.y()/2. -offset.y() +cut_extra.y(),  // y  
                          +si_dim.z()/2. -pcb_dim.z()/2. +offset.z() -cut_extra.z()); // z

    // Physical volumes

    name  = "TistarSiLayerPV";
    G4Box * si_PV = new G4Box(name,si_dim.x()/2.,si_dim.y()/2.,si_dim.z()/2.);

    G4Box * pcb_precut_PV = NULL;
    G4Box * pcb_cutting_PV = NULL;
    G4SubtractionSolid * pcb_PV = NULL;
    if(fPCBDimensionsSet) {
        name = "TistarPCBLayerPreCutPV";
        pcb_precut_PV = new G4Box(name,pcb_dim.x()/2.,pcb_dim.y()/2.,pcb_dim.z()/2.);    
        name = "TistarPCBLayerCuttingPV";
        pcb_cutting_PV = new G4Box(name,pcb_dim.x()/2.+cut_extra.x(),si_dim.y()/2.+cut_extra.y(),si_dim.z()/2.+cut_extra.z());
        name = "TistarPCBLayerPV";
        pcb_PV = new G4SubtractionSolid(name,pcb_precut_PV,pcb_cutting_PV,rotate,move);
    }

    // Logical volumes
    if(fLogicalSiLayer == NULL) {
        name = "TistarSiLayerLV";
        fLogicalSiLayer = new G4LogicalVolume(si_PV,si_material,name,0,0,0);
        fLogicalSiLayer->SetVisAttributes(si_vis_att);
    }

    if(fPCBDimensionsSet) {
        if(fLogicalPCBLayer == NULL) {
            name = "TistarPCBLayerLV";
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
            move = G4ThreeVector( 0.,                                            // x
                                  -pcb_dim.y()/2. + si_dim.y()/2. + offset.y(),  // y
                                  +pcb_dim.z()/2. - si_dim.z()/2. - offset.z()); // z
            fAssemblyLayer->AddPlacedVolume(fLogicalPCBLayer,move,rotate);
        }
    }

    return 1;
}

G4int DetectionSystemTistar::Add2StripLayer(G4double dist_from_beam, G4bool si_centered, G4LogicalVolume* expHallLog) 
{
    G4ThreeVector move;
    G4ThreeVector rotate;

    G4String name;

    TistarSettings * sett = TistarSettings::Get();

    // first strip
    move = G4ThreeVector(   +dist_from_beam,        // x
                            0.,                     // y
                            0.);                    // z
    if(fPCBDimensionsSet && !si_centered) move += G4ThreeVector(0., -fSiOffsetInPCB.y() -fSiDimensions.y()/2. +fPCBDimensions.y()/2., 0.);
    rotate = G4ThreeVector( 0.,                     // x-rotation
                            0.,                     // y-rotation
                            0.);                  // z-rotation
    name = "TistarSiLayerLV_" + std::to_string(fDetectorNumber);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);
    std::cout << "fDetectorNumber = " << fDetectorNumber << " ";
    std::cout << "layerN = " << fDetectorNumber/10%10-1 << std::endl;
    sett->SetLayerDimensionVector(fDetectorNumber/10%10-1, 0, g4_to_root(fSiDimensions));
    sett->SetLayerPositionVector(fDetectorNumber/10%10-1, 0, g4_to_root(move));

    // second strip
    move = G4ThreeVector(   -dist_from_beam,        // x
                            0.,                     // y
                            0.);                    // z
    if(fPCBDimensionsSet && !si_centered) move += G4ThreeVector(0., +fSiOffsetInPCB.y() +fSiDimensions.y()/2. -fPCBDimensions.y()/2., 0.);
    rotate = G4ThreeVector( 180,                      // x-rotation
                            0.,                  // y-rotation
                            0.);                  // z-rotation
    name = "TistarSiLayerLV_" + std::to_string(fDetectorNumber+1);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);
    std::cout << "fDetectorNumber = " << fDetectorNumber << " ";
    std::cout << "layerN = " << fDetectorNumber/10%10-1 << std::endl;
    sett->SetLayerDimensionVector(fDetectorNumber/10%10-1, 1, g4_to_root(fSiDimensions));
    sett->SetLayerPositionVector(fDetectorNumber/10%10-1, 1, g4_to_root(move));

    return 1;
}

G4int DetectionSystemTistar::Add4StripLayer(G4double dist_from_beam, G4double gap_z, G4LogicalVolume* expHallLog) 
{
    G4ThreeVector move;
    G4ThreeVector rotate;

    G4String name;
    
    TistarSettings * sett = TistarSettings::Get();

    //fPositionOffset = G4ThreeVector(0., 10.0*mm*(-5.0*mm -100.0/2.*mm +140.0/2.*mm)/(30.0*mm), 0.);

    // first forward strip
    move = G4ThreeVector(   +dist_from_beam,                // x
                            0.,                             // y
                            +fSiDimensions.z()/2. +gap_z/2.);  // z
    move += fPositionOffset;
    rotate = G4ThreeVector( 0.,                     // x-rotation
                            0.,                     // y-rotation
                            0.);                  // z-rotation
    name = "TistarSiLayerLV_" + std::to_string(fDetectorNumber);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);
    std::cout << "fDetectorNumber = " << fDetectorNumber << " ";
    std::cout << "layerN = " << fDetectorNumber/10%10-1 << std::endl;
    sett->SetLayerDimensionVector(fDetectorNumber/10%10-1, 0, g4_to_root(fSiDimensions));
    sett->SetLayerPositionVector(fDetectorNumber/10%10-1, 0, g4_to_root(move));

    // second forward strip
    move = G4ThreeVector(   -dist_from_beam,                // x
                            0.,                             // y
                            +fSiDimensions.z()/2. +gap_z/2.);  // z
    move -= fPositionOffset;
    rotate = G4ThreeVector( 0.,                     // x-rotation
                            0.,                     // y-rotation
                            180.);                  // z-rotation
    name = "TistarSiLayerLV_" + std::to_string(fDetectorNumber+1);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);
    std::cout << "fDetectorNumber = " << fDetectorNumber << " ";
    std::cout << "layerN = " << fDetectorNumber/10%10-1 << std::endl;
    sett->SetLayerDimensionVector(fDetectorNumber/10%10-1, 1, g4_to_root(fSiDimensions));
    sett->SetLayerPositionVector(fDetectorNumber/10%10-1, 1, g4_to_root(move));
    
    // first backward strip
    move = G4ThreeVector(   +dist_from_beam,                // x
                            0.,                             // y
                            -fSiDimensions.z()/2. -gap_z/2.);  // z
    move += fPositionOffset;
    rotate = G4ThreeVector( 180,                      // x-rotation
                            0.,                  // y-rotation
                            180.);                  // z-rotation
    name = "TistarSiLayerLV_" + std::to_string(fDetectorNumber+2);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);
    std::cout << "fDetectorNumber = " << fDetectorNumber << " ";
    std::cout << "layerN = " << fDetectorNumber/10%10-1 << std::endl;
    sett->SetLayerDimensionVector(fDetectorNumber/10%10-1, 2, g4_to_root(fSiDimensions));
    sett->SetLayerPositionVector(fDetectorNumber/10%10-1, 2, g4_to_root(move));

    // second backward strip
    move = G4ThreeVector(   -dist_from_beam,                // x
                            0.,                             // y
                            -fSiDimensions.z()/2. -gap_z/2.);  // z
    move -= fPositionOffset;
    rotate = G4ThreeVector( 180,                      // x-rotation
                            0.,                  // y-rotation
                            0.);                  // z-rotation
    name = "TistarSiLayerLV_" + std::to_string(fDetectorNumber+3);
    fLogicalSiLayer->SetName(name);
    PlaceDetector(move, rotate, expHallLog);
    std::cout << "fDetectorNumber = " << fDetectorNumber << " ";
    std::cout << "layerN = " << fDetectorNumber/10%10-1 << std::endl;
    sett->SetLayerDimensionVector(fDetectorNumber/10%10-1, 3, g4_to_root(fSiDimensions));
    sett->SetLayerPositionVector(fDetectorNumber/10%10-1, 3, g4_to_root(move));


    return 1;
}

// ---------------------------- Gas target ----------------------------

G4int DetectionSystemTistar::AddGasTarget(G4LogicalVolume* expHallLog)
{
    G4ThreeVector move;
    G4RotationMatrix * rotate = NULL;

    G4String target_material_name = "TistarTarget_" + TistarSettings::Get()->GetTargetMaterialName();
    G4Material * target_material = G4Material::GetMaterial(target_material_name);
    G4Material * mylar_material = G4Material::GetMaterial(fGasTargetMylarMaterialName);    
    G4Material * be_material = G4Material::GetMaterial(fGasTargetBeWindowMaterialName);    

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
        fLogicalGasTarget = new G4LogicalVolume(gas_target,target_material,"Gas_target_LV",0,0,0);
        fLogicalGasTarget->SetVisAttributes(gas_vis_att);
    }
    // Placement$a
    move = G4ThreeVector(0.,0.,0.);
    rotate = new G4RotationMatrix;
    G4VPhysicalVolume * gas_target_PV = new G4PVPlacement(rotate, move, fLogicalGasTarget, "Gas_target_PV", expHallLog, 0, 0, 0);
    
    // calculate the area density/target thickness and material density based on the simulated target parameters
    G4double outerRadius = TistarSettings::Get()->GetTargetDiameter()/2.;
    G4double halfThickness = TistarSettings::Get()->GetTargetPhysicalLength()/2;
    G4double targetThickness = fLogicalGasTarget->GetMass()/(TMath::Pi()*TMath::Power(outerRadius,2));
    G4double targetMatDensity = fLogicalGasTarget->GetMass()/(TMath::Pi()*TMath::Power(outerRadius,2)*2.*halfThickness);
    
    TistarSettings::Get()->SetTargetThickness(targetThickness/(CLHEP::mg/CLHEP::cm2));
    TistarSettings::Get()->SetTargetMaterialDensity(targetMatDensity/(CLHEP::g/CLHEP::cm3));

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

    std::cout<<"Built gas target with material \""<<target_material->GetName()<<"\", pressure "<<target_material->GetPressure()/bar*1000.<<" mbar, length "<<2.*halfThickness/cm<<" cm, "<<std::endl<<"material density "<<targetMatDensity/(g/cm3)<<" g/cm3, and area density "<<targetThickness/(mg/cm2)<<" mg/cm2"<<" calculated via target mass "<<fLogicalGasTarget->GetMass()<<std::endl;
    
    std::cout<<"Pressure = "<<target_material->GetPressure()/bar<<" bar"<<std::endl;
    std::cout<<"Temperature = "<<target_material->GetTemperature()/kelvin<<" K"<<std::endl;
    std::cout<<"Mass of molecule = "<<target_material->GetMassOfMolecule()/g<<" g"<<std::endl;
    std::cout<<"Mass of 1 mole of molecules = "<<target_material->GetMassOfMolecule()*CLHEP::Avogadro/g<<" g/mole"<<std::endl;
    std::cout<<"Density calculated via ideal gas law = "
             << ( (target_material->GetPressure()) * (target_material->GetMassOfMolecule()*CLHEP::Avogadro) / (8.314*(joule/kelvin/mole)) / (target_material->GetTemperature()) ) / (g/cm3) << " g/cm3" << std::endl;

    return 1;
}

//----------------------------- Vacuum Chamber -------------------------------

void DetectionSystemTistar::SetVacuumChamberShape(G4String shape) 
{
    if(shape=="box" || shape=="cylinder" || shape=="sphere") {
        fVacuumChamberShape = shape;
        TistarSettings::Get()->SetVacuumChamberType(shape);
    } else {
        G4cout << " ---> TI-STAR vacuum chamber shape \"" << shape << "\" unknown!" << G4endl;
    }
}

void DetectionSystemTistar::SetVacuumChamberMaterialName(G4String material) { 
    fVacuumChamberMaterialName = material; 
    TistarSettings::Get()->SetVacuumChamberGas(material);
}


G4int DetectionSystemTistar::AddVacuumChamber(G4LogicalVolume* expHallLog, G4LogicalVolume *& vacuumChamberGasLog) 
{
    TistarSettings::Get()->IncludeVacuumChamber(1);

    G4ThreeVector move;
    G4RotationMatrix * rotate = NULL;

    // Make the target materials
    G4Material * vacuum_material = G4Material::GetMaterial(fVacuumChamberMaterialName);
    G4Material * exterior_material = G4Material::GetMaterial(fVacuumChamberExteriorMaterialName);
    //TistarSettings::Get()->SetVacuumChamberGasPressure(vacuum_material->GetPressure());
    TistarSettings::Get()->SetVacuumChamberGasPressure(0.1*CLHEP::bar);

    // Set up colours and other vis. attributes
    G4VisAttributes * vacuum_vis_att = new G4VisAttributes(G4Colour::Red());
    vacuum_vis_att->SetForceWireframe(true);
    vacuum_vis_att->SetVisibility(true);

    G4VisAttributes * exterior_vis_att = new G4VisAttributes(G4Colour::Grey());
    exterior_vis_att->SetVisibility(true);
    
    // Build the object solid volume
    G4VSolid * vacuum_chamber_gas_SV = NULL;
    G4VSolid * vacuum_chamber_exterior_cutting_SV = NULL;
    G4VSolid * vacuum_chamber_exterior_SV = NULL;

    if(fVacuumChamberShape == "box") {
        // the gas inside the vacuum chamber
        vacuum_chamber_gas_SV = new G4Box("vacuum_chamber_gas_solid", 
                                      fVacuumChamberBoxDimensions.x()/2.,  // x
                                      fVacuumChamberBoxDimensions.y()/2.,  // y
                                      fVacuumChamberBoxDimensions.z()/2.); // z
        // the exterior of the vacuum chamber
        vacuum_chamber_exterior_SV = new G4Box("vacuum_chamber_exterior_precut_solid", 
                                      fVacuumChamberBoxDimensions.x()/2.+fVacuumChamberExteriorThickness,  // x
                                      fVacuumChamberBoxDimensions.y()/2.+fVacuumChamberExteriorThickness,  // y
                                      fVacuumChamberBoxDimensions.z()/2.+fVacuumChamberExteriorThickness); // z
        // hollowing out the vacuum chamber exterior
        vacuum_chamber_exterior_SV = new G4SubtractionSolid("vacuum_chamber_exterior_postcut_solid", vacuum_chamber_exterior_SV, vacuum_chamber_gas_SV); // cutting the hollow box
        // cylinder used to cut the beam exit hole in the chamber exterior
        vacuum_chamber_exterior_cutting_SV = new G4Tubs("vacuum_chamber_exterior_cutting_solid", 
                                        0., // inner radius
                                        fVacuumChamberBeamHoleRadius, // outer radius 
                                        1.*m, // z distance
                                        0, // phi start
                                        2.*M_PI); // phi end
        // cutting the beam hole in the vacuum chamber exterior
        vacuum_chamber_exterior_SV = new G4SubtractionSolid("vacuum_chamber_exterior_solid", vacuum_chamber_exterior_SV, vacuum_chamber_exterior_cutting_SV); 
    } 
    else if(fVacuumChamberShape == "cylinder") {
        // the gas inside the vacuum chamber
        vacuum_chamber_gas_SV = new G4Tubs("vacuum_chamber_gas_solid", 
                                        0., // inner radius  
                                        fVacuumChamberCylinderRadius, // outer radius 
                                        fVacuumChamberCylinderZ/2., // z distance  
                                        0, // phi start
                                        2.*M_PI); // phi end
        // the exterior of the vacuum chamber
        vacuum_chamber_exterior_SV = new G4Tubs("vacuum_chamber_exterior_precut_solid", 
                                        0., // inner radius
                                        fVacuumChamberCylinderRadius+fVacuumChamberExteriorThickness, // outer radius 
                                        fVacuumChamberCylinderZ/2.+fVacuumChamberExteriorThickness, // z distance
                                        0, // phi start
                                        2.*M_PI); // phi end
        // hollowing out the vacuum chamber exterior
        vacuum_chamber_exterior_SV = new G4SubtractionSolid("vacuum_chamber_postcut_exterior_solid", vacuum_chamber_exterior_SV, vacuum_chamber_gas_SV); // cutting the hollow box
        // cylinder used to cut the beam exit hole in the chamber exterior
        vacuum_chamber_exterior_cutting_SV = new G4Tubs("vacuum_chamber_exterior_cutting_solid", 
                                        0., // inner radius
                                        fVacuumChamberBeamHoleRadius, // outer radius 
                                        1.*m, // z distance
                                        0, // phi start
                                        2.*M_PI); // phi end
        // cutting the beam hole in the vacuum chamber exterior
        vacuum_chamber_exterior_SV = new G4SubtractionSolid("vacuum_chamber_exterior_solid", vacuum_chamber_exterior_SV, vacuum_chamber_exterior_cutting_SV); 
    } 
    else if(fVacuumChamberShape == "sphere") {
        // the gas inside the vacuum chamber
        vacuum_chamber_gas_SV = new G4Sphere("vacuum_chamber_gas_solid",
                                         0., // inner radius
                                         fVacuumChamberSphereRadius, // outer radius
                                         0., // phi start
                                         2.*M_PI, // phi end
                                         0., // theta start
                                         M_PI); // theta end
        // the exterior of the vacuum chamber
        vacuum_chamber_exterior_SV = new G4Sphere("vacuum_chamber_exterior_solid",
                                         fVacuumChamberSphereRadius, // inner radius
                                         fVacuumChamberSphereRadius+fVacuumChamberExteriorThickness, // outer radius
                                         0., // phi start
                                         2.*M_PI, // phi end
                                         0., // theta start
                                         M_PI); // theta end
        // cylinder used to cut the beam exit hole in the chamber exterior
        vacuum_chamber_exterior_cutting_SV = new G4Tubs("vacuum_chamber_exterior_cutting_solid", 
                                        0., // inner radius
                                        fVacuumChamberBeamHoleRadius, // outer radius 
                                        1.*m, // z distance
                                        0, // phi start
                                        2.*M_PI); // phi end
        // cutting the beam hole in the vacuum chamber exterior
        vacuum_chamber_exterior_SV = new G4SubtractionSolid("vacuum_chamber_exterior_solid", vacuum_chamber_exterior_SV, vacuum_chamber_exterior_cutting_SV); 
    }
    else {
        G4cout << " ---> Unknown vacuum chamber shape \"" << fVacuumChamberShape << "\", cannot build!" << G4endl;
        return 2;
    }
    
    // Logical volumes
    vacuumChamberGasLog = new G4LogicalVolume(vacuum_chamber_gas_SV, vacuum_material, "vacuum_chamber_gas_LV", 0, 0, 0);
    vacuumChamberGasLog->SetVisAttributes(vacuum_vis_att);

    G4LogicalVolume * vacuumChamberExteriorLog = new G4LogicalVolume(vacuum_chamber_exterior_SV, exterior_material, "vacuum_chamber_exterior_LV", 0, 0, 0);
    vacuumChamberExteriorLog->SetVisAttributes(exterior_vis_att);

    // Placement
    move = G4ThreeVector(0.,0.,0.);
    rotate = new G4RotationMatrix;
    G4VPhysicalVolume * vacuum_chamber_PV = new G4PVPlacement(rotate, move, vacuumChamberGasLog, "vacuum_chamber_PV", expHallLog, 0, 0, 0);    
    G4VPhysicalVolume * exterior_chamber_PV = new G4PVPlacement(rotate, move, vacuumChamberExteriorLog, "vacuum_chamber_exterior_PV", expHallLog, 0, 0, 0);

    return 1;
}

std::vector<ParticleMC>* DetectionSystemTistar::GetParticleMCvector() {
    std::vector<ParticleMC>* particleMCvector = new std::vector<ParticleMC>;
    //particleMCvector->push_back(*fBarrelErestSingleSensitiveDetector->GetParticleMC()); // THIS IS WHERE A LOT OF WORK NEEDS TO BE DONE...
    return particleMCvector;
}
