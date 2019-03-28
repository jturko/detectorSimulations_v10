/************************************************************************
 * \file MiniBallTrackInformation.hh
 *
 * \class MiniBallTrackInformation
 * \brief Storing of origin iformation of a track
 *
 * \author v.bildstein@mpi-hd.mpg.de
 *
 ************************************************************************/

#include "MiniBallTrackInformation.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

G4Allocator<MiniBallTrackInformation> aTrackInformationAllocator;

MiniBallTrackInformation::MiniBallTrackInformation()
{
    originalTraceID = 0;
    originalPosition = G4ThreeVector(0.,0.,0.);
    originalEnergy = 0.;
    originalTime = 0.;
    creatorProcessName="Init";
}

MiniBallTrackInformation::MiniBallTrackInformation(const G4Track* aTrack)
{
    originalTraceID = aTrack->GetTrackID();
    originalPosition = aTrack->GetPosition();
    originalEnergy=aTrack->GetKineticEnergy();
    originalTime = aTrack->GetGlobalTime();
    
    const G4VProcess* creatorProcess = aTrack->GetCreatorProcess();
    if(creatorProcess)
      creatorProcessName=creatorProcess->GetProcessName();
    else
      creatorProcessName="Init";
}

MiniBallTrackInformation::MiniBallTrackInformation(const MiniBallTrackInformation* aTrackInfo)
{
    originalTraceID = aTrackInfo->originalTraceID;
    originalPosition = aTrackInfo->originalPosition;
    originalEnergy = aTrackInfo->originalEnergy;
    originalTime = aTrackInfo->originalTime;
    creatorProcessName=aTrackInfo->creatorProcessName;
}

MiniBallTrackInformation::~MiniBallTrackInformation(){;}

void MiniBallTrackInformation::Print() const
{
    G4cout 
      << "Trace Nb "<<originalTraceID<< ": "
      << "CreatorProcess " << creatorProcessName 
      << " at " << originalPosition
      << "Trace Nb "<<originalTraceID<<G4endl;
}
