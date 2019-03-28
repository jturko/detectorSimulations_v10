/*
 * TRexBeam.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 *
 * Modified 2017/06/15 trockman
 * Do not call DefineNuclei until physics list is instantiated
 * 
 * Modified 2017/06/15 dhymers
 * Moved calls depending on DefineNuclei to occur afterwards
 */

#include "TRexBeam.hh"
#include "TRexSettings.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4UnitsTable.hh"

TRexBeam::TRexBeam() :
	fGammaTheta(new std::vector<G4double>(0)), fGammaPhi(new std::vector<G4double>(0)), fGammaEnergy(new std::vector<G4double>(0)),
	fGammaLab(new std::vector<G4LorentzVector>(0)) {
		// define guns
		fParticleGunEjectile = new G4ParticleGun(1);
		fParticleGunRecoil = new G4ParticleGun(1);
		fParticleGunGamma = new G4ParticleGun(1);

		fBeamEnergy = TRexSettings::Get()->GetBeamEnergy();
		fBeamWidth = TRexSettings::Get()->GetBeamWidth();
		fReactionEnergy = 0.;

		// define nuclei
		//DefineNuclei();

		// define reaction kinematics and energy loss calculations
		//fTargetMaterial = GetTargetMaterial();
		//std::cout << "TargetMaterialName for energy loss calculation in the target = " << fTargetMaterial->Name() << std::endl;
		//fKinematics = new Kinematic(&fProjectile, fTargetMaterial, TRexSettings::Get()->GetTargetThickness()/(CLHEP::mg/CLHEP::cm2));

		// energy loss in the targe
		//fEnergyVsTargetDepth = *(fKinematics->EnergyVsThickness(fBeamEnergy / CLHEP::MeV, TRexSettings::Get()->GetTargetThickness() / 1000 / (CLHEP::mg/CLHEP::cm2)));

		//	TFile bla("bla.root", "recreate");
		//	bla.cd();
		//	fEnergyVsTargetDepth.Write();
		//	bla.Close();

		// set minimal thetaCM
		fThetaCM_min = TRexSettings::Get()->GetThetaCmMin();

		fEbeamCmHist = nullptr;
		
	}

TRexBeam::~TRexBeam() {
	// TODO Auto-generated destructor stub
}

void TRexBeam::ShootReactionPosition() {
	//select random x and y position on a disk with diameter beamWidth
	/*do {
	  fReactionX = CLHEP::RandFlat::shoot(-fBeamWidth / 2., fBeamWidth / 2.) * CLHEP::mm;
	  fReactionY = CLHEP::RandFlat::shoot(-fBeamWidth / 2., fBeamWidth / 2.) * CLHEP::mm;
	  } while(sqrt(pow(fReactionX,2)+pow(fReactionY,2)) > fBeamWidth / 2.); original commented out by Leila because X/Y was not a flat distribution (gauss)*/

	fReactionX = CLHEP::RandFlat::shoot(-fBeamWidth / 2., fBeamWidth / 2.) * CLHEP::mm;
	fReactionY = CLHEP::RandFlat::shoot(-fBeamWidth / 2., fBeamWidth / 2.) * CLHEP::mm;

	// choose z according to a flat distribution in the target
	//fReactionZ = CLHEP::RandFlat::shoot(-TRexSettings::Get()->GetTargetThickness() / (2. * TRexSettings::Get()->GetTargetMaterialDensity()) / CLHEP::um,
	//TRexSettings::Get()->GetTargetThickness() / (2. * TRexSettings::Get()->GetTargetMaterialDensity()) / CLHEP::um) * CLHEP::um;
	//fReactionZ = CLHEP::RandFlat::shoot(-0.5, 0.5) * CLHEP::mm;
	fReactionZ = CLHEP::RandFlat::shoot(-TRexSettings::Get()->GetTargetPhysicalLength()/(2*CLHEP::um), TRexSettings::Get()->GetTargetPhysicalLength()/(2*CLHEP::um))*CLHEP::um;
	// units: although the target length is given as cm in the setting file but fReactionZ is in mm!
	
}

void TRexBeam::DefineNuclei() {
	fProjectileZ = TRexSettings::Get()->GetProjectileZ();
	fProjectileA = TRexSettings::Get()->GetProjectileA();
	fTargetZ = TRexSettings::Get()->GetTargetZ();
	fTargetA = TRexSettings::Get()->GetTargetA();
	fEjectileZ = TRexSettings::Get()->GetEjectileZ();
	fEjectileA = TRexSettings::Get()->GetEjectileA();
	fRecoilZ = TRexSettings::Get()->GetRecoilZ();
	fRecoilA = TRexSettings::Get()->GetRecoilA();

	// masses
	fProjectileRestMass = ParticleDefinition(fProjectileZ, fProjectileA - fProjectileZ, 0)->GetPDGMass();
	fTargetRestMass = ParticleDefinition(fTargetZ, fTargetA - fTargetZ, 0)->GetPDGMass();
	fEjectileRestMass = ParticleDefinition(fEjectileZ, fEjectileA - fEjectileZ, 0)->GetPDGMass();
	fRecoilRestMass = ParticleDefinition(fRecoilZ, fRecoilA - fRecoilZ, 0)->GetPDGMass();

	// define isotopes
	std::cout<<"Reading isotopes from '"<<TRexSettings::Get()->GetMassFile()<<"' ... ";
	fIsotopeTable = new Isotopes(TRexSettings::Get()->GetMassFile().c_str());
	if(fIsotopeTable->NumberOfIsotopes() == 0) {
		std::cout<<"failed to read mass file!"<<std::endl;
		exit(1);
	}
	std::cout<<"read "<<fIsotopeTable->NumberOfIsotopes()<<" isotopes"<<std::endl;

	fProjectile = *(fIsotopeTable->Search((char*)TRexSettings::Get()->GetProjectileName().c_str()));
	fTarget = *(fIsotopeTable->Search((char*)TRexSettings::Get()->GetTargetName().c_str()));
	//fTarget = *(IsotopeTable->Search(fTargetZ, fTargetA - fTargetZ));
	fEjectile = *(fIsotopeTable->Search((char*)TRexSettings::Get()->GetEjectileName().c_str()));
	fRecoil = *(fIsotopeTable->Search((char*)TRexSettings::Get()->GetRecoilName().c_str()));

	std::cout << "Shooting the projectile " << fProjectile.A() << fProjectile.Name() << " with (Z,A) = (" << fProjectileZ << "," <<  fProjectileA
		<< ") on the target " << fTarget.A() << fTarget.Name() << " with (Z,A) = (" << fTargetZ << "," << fTargetA << ") => ejectile "
		<< fEjectile.A() << fEjectile.Name() << " with (Z,A) = (" << fEjectileZ << "," << fEjectileA  << ") with recoil "
		<< fRecoil.A() << fRecoil.Name() << " with (Z,A) = (" << fRecoilZ << "," << fRecoilA << ")." << std::endl;

	// check settings
	if(fProjectile.Z() != fProjectileZ || fProjectile.A() != fProjectileA ||
			fTarget.Z() != fTargetZ || fTarget.A() != fTargetA ||
			fEjectile.Z() != fEjectileZ || fEjectile.A() != fEjectileA ||
			fRecoil.Z() != fRecoilZ || fRecoil.A() != fRecoilA) {
		std::cerr << "Given particle names do not match to the given charge and mass numbers!" << std::endl;
		exit(1);
	}
}

Material* TRexBeam::GetTargetMaterial() {
	Material* TargetMaterial;

	//PE and MY are implemented as materials, everything else should be the name of the element
	if(((G4String)TRexSettings::Get()->GetTargetMaterialName()).contains("PE") || ((G4String)TRexSettings::Get()->GetTargetMaterialName()).contains("MY")) {
		TargetMaterial = new Material((char*)TRexSettings::Get()->GetTargetMaterialName().c_str());
	} else {
		//if target material name is the same as the name of the scattering target build set the material to only this element
		if(TRexSettings::Get()->GetTargetMaterialName() == TRexSettings::Get()->GetTargetName() || TRexSettings::Get()->GetTargetMaterialName() == "dummy" ||
				TRexSettings::Get()->GetTargetMaterialName() == "SolidDeuterium") {       // added bei Leila 
			TargetMaterial = new Material((char*)TRexSettings::Get()->GetTargetName().c_str(),false);
		} else {
			std::cout<<"'"<<TRexSettings::Get()->GetTargetMaterialName()<<"' != '"<<TRexSettings::Get()->GetTargetName()<<"'"<<std::endl;
			char* ElementNames[] = {(char*)TRexSettings::Get()->GetTargetMaterialName().c_str(), (char*)TRexSettings::Get()->GetTargetName().c_str()};

			std::string strCarrierA = TRexSettings::Get()->GetTargetMaterialName();
			strCarrierA.erase(strCarrierA.find_first_not_of("0123456789"));
			int CarrierA = atoi(strCarrierA.c_str());

			std::string strTargetA = TRexSettings::Get()->GetTargetName();
			strTargetA.erase(strTargetA.find_first_not_of("0123456789"));
			int TargetA = atoi(strTargetA.c_str());

			G4double TargetRatio = TargetA * TRexSettings::Get()->GetTargetAtomicRatio() / (TargetA * TRexSettings::Get()->GetTargetAtomicRatio() + CarrierA);
			std::cout<<"TargetRatio = "<<TargetRatio<<" ("<<TargetA<<"*"<<TRexSettings::Get()->GetTargetAtomicRatio()<<"/("<<TargetA<<"*"<<TRexSettings::Get()->GetTargetAtomicRatio()<<"+"<<CarrierA<<"))"<<std::endl;

			double ElementRatios[] = {1-TargetRatio,TargetRatio};
			std::cout << "Element 0: " << ElementNames[0] << " with ratio " << ElementRatios[0] << std::endl;
			std::cout << "Element 1: " << ElementNames[1] << " with ratio " << ElementRatios[1] << std::endl;
			TargetMaterial = new Material(2,ElementNames,ElementRatios,false);

			//TargetMaterial = new Material((char*)TRexSettings::Get()->GetTargetMaterialName().c_str(),false);
		}
	}

	return TargetMaterial;
}

/*void TRexBeam::FillCrossSectionGraph() {
	std::ifstream file(TRexSettings::Get()->GetCrossSectionFile().c_str());

	if(file.bad()) {
		std::cerr << "Unable to open cross sectoin file" << TRexSettings::Get()->GetCrossSectionFile() << "!\nexiting ... \n";
		exit(2);
	} else {
		std::cout << "\nReading cross section file " << TRexSettings::Get()->GetCrossSectionFile() << " ... \n"<< std::endl;
	}


	// number of energies = number of lines

	//file.ignore(1000, '\n'); // ignore the first line

	file >> fNbOfBeamEnergyInCm;

	std::cout << "fNbOfBeamEnergyInCm = " << fNbOfBeamEnergyInCm << std::endl;

	// resize the vectors
	fEbeamCm.resize(fNbOfBeamEnergyInCm);
	fsigmaForEbeamCm.resize(fNbOfBeamEnergyInCm);

	// loop over all lines
	for(int i = 0; i<fNbOfBeamEnergyInCm; i++) {				
		file >>fEbeamCm[i] >>fsigmaForEbeamCm[i];
		//std::cout << "counts i: "<<i<<"	Ebeam: "<<fEbeamCm[i]<<"	fsigmaForEbeamCm: " << fsigmaForEbeamCm[i] << std::endl;				       
	} 

	file.close();

	fGraphCrossSection.push_back(TGraph(fEbeamCm.size(),&fEbeamCm[0],&fsigmaForEbeamCm[0]));
	fGrp = new TGraph(fEbeamCm.size(),&fEbeamCm[0],&fsigmaForEbeamCm[0]);
	fGrp->Draw("AL*");
	fGrp->SetMinimum(1.0e-10);
	//fGrp->SetMaximum(1.);
	fSigmaVsEbeamCmMax = fGrp->GetMaximum();
	fNumberOfPointsGraph = fGrp->GetN();
	fYaxs = fGrp->GetY();
	int locmax = TMath::LocMax(fNumberOfPointsGraph,fYaxs);
	fSigmaVsEbeamCmMax = fYaxs[locmax];		

	// *************************** converting the TGraph into a histogram ************************************

	double ebeamCmMin, ebeamCmMax;
	double sigmaEbeamCmMin, sigmaEbeamCmMax;
	double sigmaEbeamCm;

	int ebeamCmnbOfBins = fNumberOfPointsGraph*100;

	fGrp->GetPoint(0, ebeamCmMin, sigmaEbeamCmMin);
	fGrp->GetPoint(fGrp->GetN()-1, ebeamCmMax, sigmaEbeamCmMax);

	double ebeamCmbinWidth = (ebeamCmMax - ebeamCmMin) / ebeamCmnbOfBins;

	// create angular distribution histogram
	fEbeamCmHist = new TH1F("fEbeamCmHist", "fEbeamCmHist",ebeamCmnbOfBins + 1, ebeamCmMin - (ebeamCmbinWidth/2.), ebeamCmMax + (ebeamCmbinWidth/2.));
	//std::cout<<"\n ebeam min: " <<ebeamCmMin<<" ebeam max: " <<ebeamCmMax<<" sigma min: " <<sigmaEbeamCmMin<<" sigma max: " <<sigmaEbeamCmMax<<"\n bin width: "<<ebeamCmbinWidth<<" ebeamCmnbOfBins: "<<ebeamCmnbOfBins<<std::endl;

	// loop over all beam energies and fill the histogram
	for(double energy = ebeamCmMin; energy < ebeamCmMax + ebeamCmbinWidth; energy += ebeamCmbinWidth) {
		sigmaEbeamCm = fGrp->Eval(energy);
		fEbeamCmHist->Fill(energy, sigmaEbeamCm);
	}

	TFile outSigmaEbeamCmFile("sigmaEbeamCm.root", "recreate");
	outSigmaEbeamCmFile.cd();
	fGrp->Write();
	fEbeamCmHist->Write();
	outSigmaEbeamCmFile.Close();
}*/


void TRexBeam::CalculateReactionEnergyInTheTarget() {

	// *********************************** Original first reactionZ then reactionEnergy*******************************

	G4double reactionPosInTarget = fReactionZ * TRexSettings::Get()->GetTargetMaterialDensity() + TRexSettings::Get()->GetTargetThickness() / 2.;
    fReactionEnergy = fEnergyVsTargetDepth.Eval(reactionPosInTarget /(CLHEP::mg/CLHEP::cm2))*CLHEP::MeV;

	//std::cout << "fReactionZ = " << fReactionZ << " ,x = " << reactionPosInTarget /(CLHEP::mg/CLHEP::cm2) << " , E(x) = " << fReactionEnergy / CLHEP::MeV << " TargetMaterialDensity: "<<TRexSettings::Get()->GetTargetMaterialDensity()/(CLHEP::mg/CLHEP::cm3)<<std::endl; 

	// *********************************** Vinzenz first reactionEnergy then reactionZ*******************************
	
	/*fRndReaction.SetSeed(0);
	double fReacProbA = fRndReaction.Rndm();
	double fReacProbB = fRndReaction.Rndm();
	
	double fSigmaTotalBarnStdr = 2.24; // barn (32029.12/180*12.57)
	
	double fReacProb = fSigmaTotalBarnStdr * 2.81865e-3 * 0.6022 / 4.; // bar * gr/cm2 (target areal density.81865 mg/cm2)
	
	//std::cout<<"\n fReacProb: "<<fReacProb<<" fEventCounter: "<<fEventCounter<<std::endl; 
	
	if((fEventCounter-1) %1000 == 0){ 
	
	fReactionEnergyCM = fEbeamCmHist->GetRandom()/1000. * CLHEP::MeV; // MeV	

	fReactionEnergy = fReactionEnergyCM*(fTargetRestMass+fProjectileRestMass)/fTargetRestMass;//MeV

	double rangeBeam = fRangeVsBeamEnergyLeila.Eval(fBeamEnergy / CLHEP::MeV);
	double rangeReaction = fRangeVsBeamEnergyLeila.Eval(fReactionEnergy / CLHEP::MeV);

	fReactionZ = (rangeBeam-rangeReaction);	
	fReactionZ = (fReactionZ * 1000. * TRexSettings::Get()->GetTargetMaterialDensity() / (CLHEP::mg/CLHEP::cm2)*10. - TRexSettings::Get()->GetTargetPhysicalLength()/2.) * CLHEP::mm;
	
	//std::cout<<"\n fReactionZ: "<<fReactionZ<<" fReactionZ/(CLHEP::mg/CLHEP::cm2): "<<fReactionZ/(CLHEP::mg/CLHEP::cm2)<<" fReactionZ*(CLHEP::mg/CLHEP::cm2): "<<fReactionZ*(CLHEP::mg/CLHEP::cm2)<<" reaction: "<<fReacProbA<<" eventno: "<<fEventCounter<<std::endl;
	
	//std::cout<<"\n TRexSettings::Get()->GetTargetMaterialDensity(): "<<TRexSettings::Get()->GetTargetMaterialDensity()<<" TRexSettings::Get()->GetTargetMaterialDensity()/(CLHEP::mg/CLHEP::cm2): "<<TRexSettings::Get()->GetTargetMaterialDensity()/(CLHEP::mg/CLHEP::cm2)<<" TRexSettings::Get()->GetTargetMaterialDensity()*(CLHEP::mg/CLHEP::cm2): "<<TRexSettings::Get()->GetTargetMaterialDensity()*(CLHEP::mg/CLHEP::cm2)<<std::endl;
	
	//std::cout<<"\n TRexSettings::Get()->GetTargetThickness(): "<<TRexSettings::Get()->GetTargetThickness()<<" TRexSettings::Get()->GetTargetThickness()/(CLHEP::mg/CLHEP::cm2): "<<TRexSettings::Get()->GetTargetThickness()/(CLHEP::mg/CLHEP::cm2)<<" TRexSettings::Get()->GetTargetThickness()*(CLHEP::mg/CLHEP::cm2): "<<TRexSettings::Get()->GetTargetThickness()*(CLHEP::mg/CLHEP::cm2)<<std::endl;	
	
    }

    else fReactionEnergyCM = -1.0;*/

}

void TRexBeam::CreateTreeBranches() {
	if(!fTree) {
		std::cout << "\n\n\nTRexBeam: Tree doesn't exist!\n\n" << std::endl;
	}
	fTree->Branch("beamEnergy", &fBeamEnergy, "beamEnergy/D");
	fTree->Branch("beamWidth", &fBeamWidth, "beamWidth/D");
	fTree->Branch("reactionEnergy", &fReactionEnergy, "reactionEnergy/D");
	fTree->Branch("reactionEnergyCM", &fReactionEnergyCM, "reactionEnergyCM/D");
	fTree->Branch("reactionX", &fReactionX, "reactionX/D");
	fTree->Branch("reactionY", &fReactionY, "reactionY/D");
	fTree->Branch("reactionZ", &fReactionZ, "reactionZ/D");
	fTree->Branch("thetaCM", &fThetaCM, "thetaCM/D");
	fTree->Branch("ejectileTheta", &fEjectileTheta, "ejectileTheta/D");
	fTree->Branch("recoilTheta", &fRecoilTheta, "recoilTheta/D");
	fTree->Branch("ejectilePhi", &fEjectilePhi, "ejectilePhi/D");
	fTree->Branch("recoilPhi", &fRecoilPhi, "recoilPhi/D");
	fTree->Branch("ejectileEnergy", &fEjectileEnergy, "ejectileEnergy/D");
	fTree->Branch("recoilEnergy", &fRecoilEnergy, "recoilEnergy/D");
	fTree->Branch("projectileZ", &fProjectileZ, "projectileZ/I");
	fTree->Branch("projectileA", &fProjectileA, "projectileA/I");
	fTree->Branch("targetZ", &fTargetZ, "targetZ/I");
	fTree->Branch("targetA", &fTargetA, "targetA/I");
	fTree->Branch("ejectileZ", &fEjectileZ, "ejectileZ/I");
	fTree->Branch("ejectileA", &fEjectileA, "ejectileA/I");
	fTree->Branch("recoilZ", &fRecoilZ, "recoilZ/I");
	fTree->Branch("recoilA", &fRecoilA, "recoilA/I");
	fTree->Branch("scatteringProbability", &fScatteringProbability, "scatteringProbability/D");
	fTree->Branch("reaction", &fReaction, "reaction/i");

	fTree->Branch("gammaTheta", &fGammaTheta);
	fTree->Branch("gammaPhi", &fGammaPhi);
	fTree->Branch("gammaEnergy", &fGammaEnergy);
}

G4ParticleDefinition* TRexBeam::ParticleDefinition(int Z, int N, double eex) {
	if(Z+N > 4) { // create ion from ion table
		return G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z, Z+N, eex);
	} else {
		if(Z == 1 && N == 0) { // proton
			return G4Proton::ProtonDefinition();
		} else if(Z == 1 && N == 1) { // deuteron
			return G4Deuteron::DeuteronDefinition();
		} else if(Z == 1 && N == 2) { // triton
			return G4Triton::TritonDefinition();
		} else if(Z == 2 && N == 1) { // 3He
			return G4He3::He3Definition();
		} else if(Z == 2 && N == 2) { // alpha
			return G4Alpha::AlphaDefinition();
		}
	}

	std::cerr << "Error in " << __PRETTY_FUNCTION__ << "shouldn't be able to reach this stage (Z = " << Z << ", N = " << N << ")" << std::endl;
	exit(1);
}

void TRexBeam::SetEjectileGun(G4Event *anEvent) {
	if(TRexSettings::Get()->SimulateEjectiles()) {
		// particle definition
		//fParticleGunEjectile->SetParticleDefinition(ParticleDefinition(fEjectileZ, fEjectileA - fEjectileZ, fReactionEnergy)); // original
		fParticleGunEjectile->SetParticleDefinition(ParticleDefinition(fEjectileZ, fEjectileA - fEjectileZ, fExcitationEnergy));

		// emission point
		fParticleGunEjectile->SetParticlePosition(G4ThreeVector(fReactionX, fReactionY, fReactionZ));

		// set energy
		fParticleGunEjectile->SetParticleEnergy(fEjectileLab.e() - fEjectileRestMass);

		// direction
		fParticleGunEjectile->SetParticleMomentumDirection(fEjectileLab.vect());

		// generate primary vertex
		fParticleGunEjectile->GeneratePrimaryVertex(anEvent);
	}

	// set variables for the tree
	fEjectileTheta = fEjectileLab.theta() / CLHEP::radian;
	fEjectilePhi = fEjectileLab.phi() / CLHEP::radian;
	fEjectileEnergy = (fEjectileLab.e() - fEjectileRestMass) / CLHEP::keV;
}

void TRexBeam::SetRecoilGun(G4Event *anEvent) {
	// particle definition
	fParticleGunRecoil->SetParticleDefinition(ParticleDefinition(fRecoilZ, fRecoilA - fRecoilZ, 0));

	// emission point
	fParticleGunRecoil->SetParticlePosition(G4ThreeVector(fReactionX, fReactionY, fReactionZ));

	// energy
	fParticleGunRecoil->SetParticleEnergy(fRecoilLab.e() - fRecoilRestMass);

	// direction
	fParticleGunRecoil->SetParticleMomentumDirection(fRecoilLab.vect());

	// generate primary vertex
	fParticleGunRecoil->GeneratePrimaryVertex(anEvent);

	// set variables for the tree
	fRecoilTheta = fRecoilLab.theta() / CLHEP::radian;
	fRecoilPhi = fRecoilLab.phi() / CLHEP::radian;
	fRecoilEnergy = (fRecoilLab.e() - fRecoilRestMass) / CLHEP::keV;
}

void TRexBeam::SetGammaGun(G4Event *anEvent) {
	// clear old event
	fGammaTheta->resize(0);
	fGammaPhi->resize(0);
	fGammaEnergy->resize(0);

	// loop over all gammas
	for(unsigned int i = 0; i < fGammaLab->size(); i++) {
		// particle definition
		fParticleGunGamma->SetParticleDefinition(G4Gamma::GammaDefinition());

		// emission point
		fParticleGunGamma->SetParticlePosition(G4ThreeVector(fReactionX, fReactionY, fReactionZ));

		// energy
		fParticleGunGamma->SetParticleEnergy((*fGammaLab)[i].e());

		// direction
		fParticleGunGamma->SetParticleMomentumDirection((*fGammaLab)[i].vect());

		// generate primary vertex
		fParticleGunGamma->GeneratePrimaryVertex(anEvent);

		// set variables for the tree
		fGammaTheta->push_back((*fGammaLab)[i].theta() / CLHEP::radian);
		fGammaPhi->push_back((*fGammaLab)[i].phi() / CLHEP::radian);
		fGammaEnergy->push_back((*fGammaLab)[i].e() / CLHEP::keV);

		//std::cout << "fGammaEnergy[" << i << "] = " << (*fGammaEnergy)[i] << std::endl;
	}
}
