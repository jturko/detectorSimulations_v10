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

DetectionSystemTISTAR::DetectionSystemTISTAR()
    // LogicalVolumes
{
}

DetectionSystemTISTAR::~DetectionSystemTISTAR() 
{
    // LogicalVolumes
}

G4int DetectionSystemTISTAR::Build() 
{
    fAssemblyTISTAR = new G4AssemblyVolume(); 

    BuildTISTAR();

    return 1;
}

G4int DetectionSystemTISTAR::PlaceDetector(G4LogicalVolume* expHallLog) {
    G4RotationMatrix * rotate = new G4RotationMatrix;
    G4ThreeVector move = G4ThreeVector(0., 0., 0.);

    fAssemblyTISTAR->MakeImprint(expHallLog, move, rotate);

    return 1;
}

G4int DetectionSystemTISTAR::BuildTISTAR() 
{
    G4ThreeVector move, direction;
    G4RotationMatrix* rotate;

    return 1;
}

