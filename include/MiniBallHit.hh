//////////////////////////////////////////////////////////////////////////////
/**  \file
 *   \brief Hit in a detector
 */ 
///////////////////////////////////////////////////////////////////////////////

#ifndef MiniBallHit_h
#define MiniBallHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include "Randomize.hh"

extern unsigned int zaehler;

/// Hit in a detector
class MiniBallHit : public G4VHit
{
public:

  MiniBallHit();
 ~MiniBallHit();
  MiniBallHit(const MiniBallHit&);
  const MiniBallHit& operator=(const MiniBallHit&);
  int operator==(const MiniBallHit&) const;
  
  inline void *operator new(size_t);
  inline void  operator delete(void *aHit);

  void Draw();
  void Print();

public:
  
  void SetTraceID   (G4int trace)      { traceID = trace; };
  void SetDetectorID(G4int det)        { detectorID = det; };
  void SetEdep      (G4double de)      { edep = de; };
  void SetPos       (G4ThreeVector xyz){ pos = xyz; };
  void SetPosInDet  (G4ThreeVector xyz){ posInDet = xyz; };
  void SetProcName  (G4String name)    { procName = name; };
  void SetSegmentID (G4int id)         { segmentID = id; };
  void SetTime      (G4double time)    { fTime = time; }//Tomer time function
  
  G4int GetTraceID()          { return traceID; };
  G4int GetDetectorID()       { return detectorID; };
  G4double GetEdep()          { return edep; };
  G4ThreeVector GetPos()      { return pos; };
  G4ThreeVector GetPosInDet() { return posInDet; };
  G4String GetProcName()      { return procName; };
  G4int GetSegmentID()        { return segmentID; };
  G4double GetTime()	      { return fTime; }//Tomer time function
  
private:
  
  G4int         traceID, detectorID;
  G4double      edep;
  G4ThreeVector pos, posInDet;
  G4String      procName;
  G4int         segmentID;
  G4double 	fTime;
};

typedef G4THitsCollection<MiniBallHit> MiniBallHitsCollection;

extern G4Allocator<MiniBallHit> MiniBallHitAllocator;

inline void* MiniBallHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) MiniBallHitAllocator.MallocSingle();
  return aHit;
}

inline void MiniBallHit::operator delete(void *aHit)
{
  MiniBallHitAllocator.FreeSingle((MiniBallHit*) aHit);
}

#endif


