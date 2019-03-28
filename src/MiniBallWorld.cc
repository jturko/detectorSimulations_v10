#include "MiniBallWorld.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "MiniBallMaterial.hh"

//include-files for visualisation
#ifdef G4VIS_USE
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#endif

MiniBallWorld::MiniBallWorld(G4double wl)
  :worldlength(wl), mother(NULL)
{;}

MiniBallWorld::~MiniBallWorld()
{;}

G4VPhysicalVolume* MiniBallWorld::Construct()
{
  G4Material* vacuum = MiniBallMaterial::Get()->GetMaterialVacuum();

  G4Box* world_sol
    = new G4Box("World_Sol",worldlength, worldlength,worldlength);

  G4LogicalVolume* world_log
    = new G4LogicalVolume(world_sol,vacuum, "World_Log"); 

#ifdef G4VIS_USE
  // Visualization attributes for world volume
  world_log->SetVisAttributes (G4VisAttributes::Invisible);
#endif

  mother = new G4PVPlacement(0, G4ThreeVector(), "World",
			     world_log, NULL, false, 0);

  return mother;
}

