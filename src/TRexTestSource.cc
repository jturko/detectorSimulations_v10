/*
 * TRexTestSource.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#include "TRexTestSource.hh"
#include "TRexSettings.hh"

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "g4root.hh"

TRexTestSource::TRexTestSource() {
	fParticleGun = new G4ParticleGun(1);
}

TRexTestSource::~TRexTestSource() {
	// TODO Auto-generated destructor stub
}

/************************************************************
 *
 * Simulates a test source
 *
 ************************************************************/
void TRexTestSource::GeneratePrimaries(G4Event *anEvent) {

	fReactionEnergy = TRexSettings::Get()->GetTestSourceEnergy();

	// shoot the alpha emission point
	ShootReactionPosition();

	fParticleGun->SetParticleDefinition(G4Proton::ProtonDefinition());
	fParticleGun->SetParticlePosition(G4ThreeVector(fReactionX, fReactionY, fReactionZ));
	fParticleGun->SetParticleEnergy(fReactionEnergy);

	// isotropic distribution
	CreateIsotropicDistribution();

	G4ThreeVector direction;
	direction.set(1,1,1);
	direction.setMag(1.);
	direction.setTheta(fThetaCM);
	direction.setPhi(fPhi);
	//direction.setTheta(45* degree);
	//direction.setPhi(2.0* degree);
	fParticleGun->SetParticleMomentumDirection(direction);

	fParticleGun->GeneratePrimaryVertex(anEvent);

	FillNtuple();
}

void TRexTestSource::ShootReactionPosition() {
	fReactionX = 0.0 * mm;
	fReactionY = 0.0 * mm;
	fReactionZ = 0.0 * mm;
}

void TRexTestSource::CreateIsotropicDistribution() {
	fThetaCM = G4RandFlat::shoot(-1., 1.);

	fThetaCM = acos(fThetaCM)*radian;

	fPhi = G4RandFlat::shoot(-M_PI,M_PI)*radian;
	//fPhi = RandFlat::shoot(0., 2.* M_PI)*radian;
	//fPhi = RandFlat::shoot(-M_PI / 2.,M_PI + M_PI / 2.)*radian;
}

void TRexTestSource::CreateNtupleBranches() {
    // w/ Griffinv10
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    fNtupleID = analysisManager->CreateNtuple("treeGen", "generator output from TRexTestSource");
    fNtupleColID[0] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionEnergy");
    fNtupleColID[1] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionX");
    fNtupleColID[2] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionY");
    fNtupleColID[3] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionZ");
    fNtupleColID[4] = analysisManager->CreateNtupleDColumn(fNtupleID, "thetaCM");
    fNtupleColID[5] = analysisManager->CreateNtupleDColumn(fNtupleID, "phi");
    analysisManager->FinishNtuple(fNtupleID);
    G4cout << "created ntuple treeGen" << G4endl;
}

void TRexTestSource::FillNtuple() {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[0], fReactionEnergy);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[1], fReactionX);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[2], fReactionY);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[3], fReactionZ);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[4], fThetaCM);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[5], fPhi);
    analysisManager->AddNtupleRow(fNtupleID);
}


