/*************************************************************************
 *
 *    hit in a segmented HPGe detector
 * 
 *************************************************************************/

#include "MiniBallHit.hh"

#include "MBglobals.hh"

#include "G4UnitsTable.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<MiniBallHit> MiniBallHitAllocator;

MiniBallHit::MiniBallHit() 
{;}

MiniBallHit::~MiniBallHit() 
{;}

MiniBallHit::MiniBallHit(const MiniBallHit& right)
  : G4VHit()
{
  traceID    = right.traceID;
  detectorID = right.detectorID;
  edep       = right.edep;
  pos        = right.pos;
  posInDet   = right.posInDet;
  procName   = right.procName;
  segmentID  = right.segmentID;
}

const MiniBallHit& MiniBallHit::operator=(const MiniBallHit& right)
{
  traceID    = right.traceID;
  detectorID = right.detectorID;
  edep       = right.edep;
  pos        = right.pos;
  posInDet   = right.posInDet;
  procName   = right.procName;
  segmentID  = right.segmentID;
  return *this;
}

int MiniBallHit::operator==(const MiniBallHit& right) const
{
  return (this==&right) ? 1 : 0;
}

void MiniBallHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void MiniBallHit::Print()
{
  G4cout << " MiniBall: " << G4endl;
  G4cout << "  traceID: " << traceID 
         << " energy deposit: " << G4BestUnit(edep, "Energy")
	 << " position: " << G4BestUnit(pos, "Length")
	 << " proc name: " << procName << G4endl;
  G4cout << "  detectorID: " << detectorID 
    	 << "segmentNb: " << segmentID  
	 << " position in detector: " << G4BestUnit(posInDet, "Length")
         << G4endl;
}


