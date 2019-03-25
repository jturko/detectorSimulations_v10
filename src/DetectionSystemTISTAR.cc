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
    fNLayers(nlayers)
{
    fSiliconMaterialName = "Silicon";
    fPCBMaterialName = "Epoxy-Resin";
   
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
    for(int i=0; i<int(fLogicalSiLayers.size()); ++i)   delete (fLogicalSiLayers[i]);
    for(int i=0; i<int(fLogicalPCBLayers.size()); ++i)  delete (fLogicalPCBLayers[i]);
    for(int i=0; i<int(fAssemblyLayers.size()); ++i)    delete (fAssemblyLayers[i]);
}

G4int DetectionSystemTISTAR::Build() 
{
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
    BuildLayer(0);

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

