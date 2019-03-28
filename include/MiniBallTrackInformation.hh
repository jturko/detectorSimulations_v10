/************************************************************************
 * \file MiniBallTrackInformation.hh
 *
 * \class MiniBallTrackInformation
 * \brief Storing of origin iformation of a track
 *
 * \author v.bildstein@mpi-hd.mpg.de
 *
 ************************************************************************/

#ifndef MiniBallTrackInformation_h
#define MiniBallInformation_h 

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class MiniBallTrackInformation : public G4VUserTrackInformation 
{
public:
  MiniBallTrackInformation();
  MiniBallTrackInformation(const G4Track* aTrack);
  MiniBallTrackInformation(const MiniBallTrackInformation* aTrackInfo);
  virtual ~MiniBallTrackInformation();
   
  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);
  inline int operator ==(const MiniBallTrackInformation& right) const
  {return (this==&right);}
  
  void Print() const;
  
private:
  G4int                 originalTraceID;
  G4ParticleDefinition* particleDefinition;
  G4ThreeVector         originalPosition;
  G4double              originalEnergy;
  G4double              originalTime;
  G4String              creatorProcessName;

  public:
  inline G4int GetOriginalTraceID() const {return originalTraceID;}
  inline G4ParticleDefinition* GetOriginalParticle() const {return particleDefinition;}
  inline G4ThreeVector GetOriginalPosition() const {return originalPosition;}
  inline G4double GetOriginalEnergy() const {return originalEnergy;}
  inline G4double GetOriginalTime() const {return originalTime;}
  inline G4String GetCreatorProcess() const {return creatorProcessName;}
};

extern G4Allocator<MiniBallTrackInformation> aTrackInformationAllocator;

inline void* MiniBallTrackInformation::operator new(size_t)
{ void* aTrackInfo;
 aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void MiniBallTrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator.FreeSingle((MiniBallTrackInformation*)aTrackInfo);}

#endif



