/*
 * TRexBeamIn.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 * 
 * modified to build projectile based on settings file Z, A
 * dhymers 2017/06/12
 */

#include "TRexBeamIn.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "TistarSettings.hh"
#include "G4ParticleGun.hh"
#include "G4Alpha.hh"
#include "Randomize.hh"

#include "TistarSettings.hh"

#include "g4root.hh"

TRexBeamIn::TRexBeamIn() {
	fParticleGun = new G4ParticleGun(1);
}

TRexBeamIn::~TRexBeamIn() {
	// TODO Auto-generated destructor stub
}

/************************************************************
 *
 * Simulates a quadruple alpha source
 *
 ************************************************************/
void TRexBeamIn::GeneratePrimaries(G4Event *anEvent) {
	
	
	ShootReactionPosition();
	fReactionEnergy = TistarSettings::Get()->GetBeamEnergy();  
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();//Creating Projectile Particle 
	G4ParticleDefinition* beam = particleTable->GetIonTable()->GetIon(
		TistarSettings::Get()->GetProjectileZ(), //Atomic Number
		TistarSettings::Get()->GetProjectileA(), //Atomic Mass
		0                               //Excitation Energy
	);
	fParticleGun->SetParticleDefinition(beam);
	fParticleGun->SetParticlePosition(G4ThreeVector(fReactionX, fReactionY, fReactionZ));
	fParticleGun->SetParticleEnergy(fReactionEnergy);
	
	// isotropic distribution, if reaction_z < 0 shoot only to negative z's (theta < 90) otherwise only to positive z's (theta > 90)
	CreateIsotropicDistribution();

	G4ThreeVector direction;
	direction.set(0,0,1);
	direction.setMag(1.);
	direction.setTheta(fThetaCM);
	direction.setPhi(fPhi);
	fParticleGun->SetParticleMomentumDirection(direction);

	fParticleGun->GeneratePrimaryVertex(anEvent);
	//std::cout << "***********fBeamEnergy***************"<<fReactionEnergy<<std::endl;
	FillNtuple();
}

void TRexBeamIn::ShootReactionPosition() {
	G4double BeamDiameter = TistarSettings::Get()->GetBeamWidth() / mm;
	fReactionZ = TistarSettings::Get()->GetGasTargetLength() / cm ;
	fReactionZ = -(fReactionZ/2.) - 0.0008;
	fReactionZ *= cm;
	
	//hard-code start position
	//fReactionZ = -1 * cm; ???

	do {
		fReactionX = CLHEP::RandFlat::shoot(-BeamDiameter / 2., BeamDiameter / 2.);
		fReactionY = CLHEP::RandFlat::shoot(-BeamDiameter / 2., BeamDiameter / 2.);
	} while(sqrt(pow(fReactionX,2) + pow(fReactionY,2)) > BeamDiameter / 2.);

	fReactionX *= mm;
	fReactionY *= mm;
}

void TRexBeamIn::CreateIsotropicDistribution() {
	//fThetaCM = CLHEP::RandFlat::shoot(0., M_PI/180.)*radian;
        fThetaCM =(0.0*M_PI/180.)*radian;
	fPhi = CLHEP::RandFlat::shoot(-M_PI,M_PI)*radian;//M_PI=180 deg
	//fPhi = CLHEP::RandFlat::shoot(0., 2.* M_PI)*radian;
	//fPhi = CLHEP::RandFlat::shoot(-M_PI / 2.,M_PI + M_PI / 2.)*radian;
}

void TRexBeamIn::CreateNtupleBranches(TTree * tree) {
    if(TistarSettings::Get()->SaveMe()) {
        fTree = tree;
        fTree->Branch("reactionEnergy", &fReactionEnergy, "reactionEnergy/D");
        fTree->Branch("reactionX", &fReactionX, "reactionX/D");
        fTree->Branch("reactionY", &fReactionY, "reactionY/D");
        fTree->Branch("reactionZ", &fReactionZ, "reactionZ/D");
        fTree->Branch("thetaCM", &fThetaCM, "thetaCM/D");
        fTree->Branch("phi", &fPhi, "phi/D");
    }
    else {
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        fNtupleID = analysisManager->CreateNtuple("treeGen", "generator output from TRexBeamIn");
        fNtupleColID[0] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionEnergy");
        fNtupleColID[1] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionX");
        fNtupleColID[2] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionY");
        fNtupleColID[3] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionZ");
        fNtupleColID[4] = analysisManager->CreateNtupleDColumn(fNtupleID, "thetaCM");
        fNtupleColID[5] = analysisManager->CreateNtupleDColumn(fNtupleID, "phi");
        analysisManager->FinishNtuple(fNtupleID);
    }
    G4cout << "created TRexBeamIn ntuple: "<<fTree->GetName()<<", "<<fTree<<G4endl;
}

void TRexBeamIn::FillNtuple() {
    if(TistarSettings::Get()->SaveMe()) {
        fTree->Fill();
    }
    else {
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[0], fReactionEnergy);
        analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[1], fReactionX);
        analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[2], fReactionY);
        analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[3], fReactionZ);
        analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[4], fThetaCM);
        analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[5], fPhi);
        analysisManager->AddNtupleRow(fNtupleID);
    }
}

