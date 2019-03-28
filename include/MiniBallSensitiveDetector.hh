///////////////////////////////////////////////////////////////////////////////
/** \file 
 *  \brief Sensitive detector for detector read-out
 **/
//////////////////////////////////////////////////////////////////////////////

#ifndef MiniBallSensitiveDetector_h
#define MiniBallSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "MiniBallHit.hh"

class G4Step;
class G4HCofThisEvent;
class MiniBallHistoManager;
class G4VPhysicalVolume;
class G4AffineTransform;


///Sensitive detector for detector read-out
class MiniBallSensitiveDetector : public G4VSensitiveDetector
{

public:
  MiniBallSensitiveDetector(G4String name);
  MiniBallSensitiveDetector(MiniBallHistoManager* HM, G4String name);
  ~MiniBallSensitiveDetector();

  void   Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void   EndOfEvent(G4HCofThisEvent*);
  
  G4int GetDetectorID() { return detectorID; };

private:
  MiniBallHitsCollection* hitsCollection;
  G4String HCname;

  MiniBallHistoManager* histoManager;

  G4int detectorID;
  G4int tracecounter;

  G4int HCID;
};

#endif

