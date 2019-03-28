/***********************************************************************
 *
 *    sensitive detector for detector read-out of each segment
 *
 ***********************************************************************/

#include "MiniBallSensitiveDetector.hh"
#include "MiniBallHistoManager.hh"
#include "MiniBallTrackInformation.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"

#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"

#include "G4VPhysicalVolume.hh"

#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

MiniBallSensitiveDetector::MiniBallSensitiveDetector(G4String name)
  :G4VSensitiveDetector(name), histoManager(0), detectorID(-1),
   tracecounter(2), HCID(-1)
{
  HCname = MB_HITS_COLLECTION_NAME;
  collectionName.insert(HCname);
}

// Initialize HitsCollection and histogramming
MiniBallSensitiveDetector::MiniBallSensitiveDetector(MiniBallHistoManager* HM,
						     G4String name)
  :G4VSensitiveDetector(name), histoManager(HM),
   detectorID(-1), tracecounter(2), HCID(-1)
{
  HCname = MB_HITS_COLLECTION_NAME;
  collectionName.insert(HCname);
  if(histoManager)
    detectorID = histoManager->AddDetector(name);
}

MiniBallSensitiveDetector::~MiniBallSensitiveDetector()
{ ;}

void MiniBallSensitiveDetector::Initialize(G4HCofThisEvent* evnt)
{
  hitsCollection = new MiniBallHitsCollection(SensitiveDetectorName, HCname); 

  if(HCID < 0)
    {
      G4String full_name = SensitiveDetectorName;
      full_name += "/" + HCname;
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(full_name);
    }

  evnt->AddHitsCollection(HCID, hitsCollection);

  if(histoManager)
    histoManager->SetMBHitsCollectionID(detectorID, HCID);
}

G4bool MiniBallSensitiveDetector::ProcessHits(G4Step* aStep,
					      G4TouchableHistory*)
{
  
  MiniBallTrackInformation* info = 
    (MiniBallTrackInformation*)(aStep->GetTrack()->GetUserInformation());

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4String process
    = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

  if(!info||process=="phot"||process=="compt"||process=="conv")
    tracecounter++;
     
  if(edep == 0.)
    return false;
 
  MiniBallHit* newHit = new MiniBallHit();
  newHit->SetEdep(edep);
  if(info)
    newHit->SetTraceID(info->GetOriginalTraceID());
  else
    newHit->SetTraceID(tracecounter);
 
  
  newHit->SetDetectorID(detectorID);
  newHit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());

  G4StepPoint* stepPoint;
  if(process=="Transportation")
    stepPoint = aStep->GetPreStepPoint();
  else
    stepPoint = aStep->GetPostStepPoint();

  G4ThreeVector position(stepPoint->GetPosition());

  newHit->SetPos(position);

  G4TouchableHandle theTouchable = stepPoint->GetTouchableHandle();

  G4ThreeVector det_position
    = theTouchable->GetHistory()->GetTopTransform().TransformPoint(position);
  newHit->SetPosInDet(det_position);

  G4double phi = det_position.phi() + 30.*deg;
  while(phi<0) phi += 2*M_PI;
  while(phi>2*M_PI) phi -= 2*M_PI;

  G4int segment_id = (G4int) (phi/(60.*deg));

  newHit->SetSegmentID(segment_id);
  
  newHit->SetProcName(process);
  
  // newHit->Print();

  hitsCollection->insert(newHit);
  
#ifdef G4VIS_USE
  // newHit->Draw();
#endif

  return true;
}

void MiniBallSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{tracecounter=1;}

	      









