//////////////////////////////////////////////////////////////////////////////
/** \file 
 *  \brief Creation of a world
 */
//////////////////////////////////////////////////////////////////////////////

#ifndef MiniBallWorld_H
#define MiniBallWorld_H 1

#include "G4VUserDetectorConstruction.hh"
#include "MBglobals.hh"

class G4VPhysicalVolume;

/// Creation of the World
class MiniBallWorld:public G4VUserDetectorConstruction
{
public:
  MiniBallWorld(G4double wl=WORLD_LENGTH);
  ~MiniBallWorld();

  G4VPhysicalVolume* Construct();

private:
  G4double worldlength;
  G4VPhysicalVolume* mother;
};
#endif
