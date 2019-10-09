/*
 * TRexAlphaSource.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#include "TRexAlphaSource.hh"

#include "TistarSettings.hh"
#include "G4ParticleGun.hh"
#include "G4Alpha.hh"
#include "G4Proton.hh"
#include "Randomize.hh"

#include "g4root.hh"

TRexAlphaSource::TRexAlphaSource() {
	fParticleGun = new G4ParticleGun(1);
}

TRexAlphaSource::~TRexAlphaSource() {
	// TODO Auto-generated destructor stub
}

/************************************************************
 *
 * Simulates a quadruple alpha source
 *
 ************************************************************/
void TRexAlphaSource::GeneratePrimaries(G4Event *anEvent) {
	//choose the emitted alpha
	double tmp = G4RandFlat::shoot(0.,3.); 

	/*if(tmp < 1) { //239Pu #####
		tmp = G4RandFlat::shoot(0., 99.82);//70.77+17.11+11.94

		if(tmp < 70.77) {
			fReactionEnergy = 5156.59*keV;
		} else if(tmp < 70.77+17.11) {
			fReactionEnergy = 5144.3*keV;
		} else {
			fReactionEnergy = 5105.5*keV;
		}
	} else if(tmp < 2) {//241Am
		tmp = G4RandFlat::shoot(0., 97.9);//84.8+13.1

		if(tmp < 84.8) {
			fReactionEnergy = 5485.56*keV;
		} else {
			fReactionEnergy = 5442.80*keV;
		}
	} else {//244Cm
		tmp = G4RandFlat::shoot(0., 100.);//76.4+23.6

		if(tmp < 76.4) {
			fReactionEnergy = 5804.77*keV;
		} else {
			fReactionEnergy = 5762.64*keV;
		}
	} ##### */
	
	fReactionEnergy = 10000.*keV; // leila
	//fReactionEnergy = G4RandFlat::shoot(0.,10000.)*keV; // leila

	// shoot the alpha emission point
	ShootReactionPosition();

	//fParticleGun->SetParticleDefinition(G4Alpha::AlphaDefinition()); 
	fParticleGun->SetParticleDefinition(G4Proton::ProtonDefinition());// leila
	fParticleGun->SetParticlePosition(G4ThreeVector(fReactionX, fReactionY, fReactionZ));
	fParticleGun->SetParticleEnergy(fReactionEnergy);

	// isotropic distribution, if reaction_z < 0 shoot only to negative z's (theta < 90) otherwise only to positive z's (theta > 90)
	CreateIsotropicDistribution();

	G4ThreeVector direction;
	direction.set(1,1,1);
	direction.setMag(1.);
	direction.setTheta(fThetaCM);
	direction.setPhi(fPhi);
	fParticleGun->SetParticleMomentumDirection(direction);

	fParticleGun->GeneratePrimaryVertex(anEvent);

	FillNtuple();
}

void TRexAlphaSource::ShootReactionPosition() {
	G4double alphaSourceDiameter = TistarSettings::Get()->GetAlphaSourceDiameter() / mm;
	G4double alphaSourceThickness = TistarSettings::Get()->GetAlphaSourceThickness();
	G4double BeamDiameter = TistarSettings::Get()->GetBeamWidth() / mm;

	do {
		/**fReactionX = G4RandFlat::shoot(-alphaSourceDiameter / 2., alphaSourceDiameter / 2.);
		fReactionY = G4RandFlat::shoot(-alphaSourceDiameter / 2., alphaSourceDiameter / 2.);
		} //while(sqrt(pow(fReactionX,2) + pow(fReactionY,2)) > alphaSourceDiameter / 2.);**/
		
		fReactionX = G4RandFlat::shoot(-BeamDiameter / 2., BeamDiameter / 2.);
		fReactionY = G4RandFlat::shoot(-BeamDiameter / 2., BeamDiameter / 2.);
		
	} while(sqrt(pow(fReactionX,2) + pow(fReactionY,2)) > BeamDiameter / 2.);

	fReactionX *= mm;
	fReactionY *= mm;

	// alpha source: only from the surface of the 'target' (i.e. the source)
	double tmp = G4RandFlat::shoot(0.,2.);

	/**if(tmp < 1.) {
		fReactionZ = -alphaSourceThickness;
	} else {
		fReactionZ = alphaSourceThickness;
	}**/
	
	//fReactionZ = TistarSettings::Get()->GetGasTargetLength() / cm ; // leila
	//fReactionZ = (fReactionZ/2.) * cm; // leila
	
	fReactionZ = G4RandFlat::shoot(-TistarSettings::Get()->GetGasTargetLength()/2.,TistarSettings::Get()->GetGasTargetLength()/2.)*mm; // leila #####
	
	//fReactionZ = alphaSourceThickness;//Leila
	
}

void TRexAlphaSource::CreateIsotropicDistribution() {
	/**if(fReactionZ < 0) {
		fThetaCM = G4RandFlat::shoot(-1., 0.);
	} else {
		fThetaCM = G4RandFlat::shoot(0., 1.);
	}
	
	fThetaCM = acos(fThetaCM)*radian;**/
	//fThetaCM =(105.0*M_PI/180.)*radian; // leila
	fThetaCM =G4RandFlat::shoot(-M_PI,M_PI)*radian;

	//fPhi = G4RandFlat::shoot(104.*M_PI/180.,106.*M_PI/180.)*radian; // leila
	//fPhi = G4RandFlat::shoot(90.*M_PI/180.)*radian; // leila
	//fPhi = G4RandFlat::shoot(-M_PI,M_PI)*radian;
	fPhi = G4RandFlat::shoot(0., 2.* M_PI)*radian;
	//fPhi = G4RandFlat::shoot(-M_PI / 2.,M_PI + M_PI / 2.)*radian;
}

void TRexAlphaSource::CreateNtupleBranches(TTree * tree) {
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
        fNtupleID = analysisManager->CreateNtuple("treeGen", "generator output from TRexAlphaSource");
        fNtupleColID[0] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionEnergy");
        fNtupleColID[1] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionX");
        fNtupleColID[2] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionY");
        fNtupleColID[3] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionZ");
        fNtupleColID[4] = analysisManager->CreateNtupleDColumn(fNtupleID, "thetaCM");
        fNtupleColID[5] = analysisManager->CreateNtupleDColumn(fNtupleID, "phi");
        analysisManager->FinishNtuple(fNtupleID);
    }     
    G4cout << "created ntuple treeGen" << G4endl;
}

void TRexAlphaSource::FillNtuple() {
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

