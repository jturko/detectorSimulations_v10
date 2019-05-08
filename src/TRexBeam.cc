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

#include "Randomize.hh"

#include "g4root.hh"

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
		//fKinematics = new Kinematic(&fProjectile, fTargetMaterial, TRexSettings::Get()->GetTargetThickness()/(mg/cm2));

		// energy loss in the targe
		//fEnergyVsTargetDepth = *(fKinematics->EnergyVsThickness(fBeamEnergy / MeV, TRexSettings::Get()->GetTargetThickness() / 1000 / (mg/cm2)));

		// set minimal thetaCM
		fThetaCM_min = TRexSettings::Get()->GetThetaCmMin();

		//fEbeamCmHist = nullptr;
        
        fSplineReactionZvsRadius = NULL;
		
	}

TRexBeam::~TRexBeam() {
	// TODO Auto-generated destructor stub
}

void TRexBeam::ShootReactionPosition() {
	// choose z according to a flat distribution in the target
    fReactionZ = G4RandFlat::shoot(-TRexSettings::Get()->GetTargetPhysicalLength()/(2*um), TRexSettings::Get()->GetTargetPhysicalLength()/(2*um))*um;
    
    // if the beam spread file has been set, use the spline to get the x-y reaction coords
	if(TRexSettings::Get()->GetSplineReactionZvsRadiusFileBool()) {
        // calculate the radius 
        double max_radius = fSplineReactionZvsRadius->Eval(fReactionZ/mm)*mm;

        double theta = 2. * M_PI * G4UniformRand();
        double radius = max_radius * sqrt( G4UniformRand() );

        fReactionX = radius * cos(theta);
        fReactionY = radius * sin(theta);
    } 
    else {
        //select random x and y position on a disk with diameter beamWidth
	    /*do {
	      fReactionX = G4RandFlat::shoot(-fBeamWidth / 2., fBeamWidth / 2.) * mm;
	      fReactionY = G4RandFlat::shoot(-fBeamWidth / 2., fBeamWidth / 2.) * mm;
	      } while(sqrt(pow(fReactionX,2)+pow(fReactionY,2)) > fBeamWidth / 2.); original commented out by Leila because X/Y was not a flat distribution (gauss)*/

	    fReactionX = G4RandFlat::shoot(-fBeamWidth / 2., fBeamWidth / 2.) * mm;
	    fReactionY = G4RandFlat::shoot(-fBeamWidth / 2., fBeamWidth / 2.) * mm;

	    // choose z according to a flat distribution in the target
	    //fReactionZ = G4RandFlat::shoot(-TRexSettings::Get()->GetTargetThickness() / (2. * TRexSettings::Get()->GetTargetMaterialDensity()) / um,
	    //TRexSettings::Get()->GetTargetThickness() / (2. * TRexSettings::Get()->GetTargetMaterialDensity()) / um) * um;
	    //fReactionZ = G4RandFlat::shoot(-0.5, 0.5) * mm;
	    //fReactionZ = G4RandFlat::shoot(-TRexSettings::Get()->GetTargetPhysicalLength()/(2*um), TRexSettings::Get()->GetTargetPhysicalLength()/(2*um))*um;
	    // units: although the target length is given as cm in the setting file but fReactionZ is in mm!
    }	

}

void TRexBeam::BuildSplineReactionZvsRadius() {
    std::ifstream file(TRexSettings::Get()->GetSplineReactionZvsRadiusFile().c_str());
    if(!TRexSettings::Get()->GetSplineReactionZvsRadiusFileBool()) {
        std::cout<<"No beam spread file set, fSplineReactionZvsRadius could not be built!\nexiting ... \n";
        exit(2);
    } else if(file.bad()) {
        std::cerr << "Unable to open beam spread distribution file" << TRexSettings::Get()->GetSplineReactionZvsRadiusFile() << "!\nexiting ... \n";
        exit(2);
    } else {
        std::cout << "\nReading beam spread distribution file " << TRexSettings::Get()->GetSplineReactionZvsRadiusFile() << " ... \n"<< std::endl;
    }    

    int nbPoints;
    file >> nbPoints;
    std::cout << "nbPoints = " << nbPoints << " ...\n" << std::endl;
    
    // build arrays 
    double * reactionZArray = new double[(const int)nbPoints];
    double * radiusArray = new double[(const int)nbPoints];

    // fill with the data
    for(int i=0; i<nbPoints; i++) {
        file >> reactionZArray[i] >> radiusArray[i];
    }

    fSplineReactionZvsRadius = new TSpline3("fSplineReactionZvsRadius", reactionZArray, radiusArray, nbPoints);

    file.close();
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
    fReactionEnergy = fEnergyVsTargetDepth.Eval(reactionPosInTarget /(mg/cm2))*MeV;

	//std::cout << "fReactionZ = " << fReactionZ << " ,x = " << reactionPosInTarget /(mg/cm2) << " , E(x) = " << fReactionEnergy / MeV << " TargetMaterialDensity: "<<TRexSettings::Get()->GetTargetMaterialDensity()/(mg/cm3)<<std::endl; 

	// *********************************** Vinzenz first reactionEnergy then reactionZ*******************************
	
	/*fRndReaction.SetSeed(0);
	//double fReacProbA = fRndReaction.Rndm(); // so we can get rid of ROOT, G4RandFlat was added
	//double fReacProbB = fRndReaction.Rndm();
    double fReacProbA = G4RandFlat::shoot(0,1);	
    double fReacProbB = G4RandFlat::shoot(0,1);	

	double fSigmaTotalBarnStdr = 2.24; // barn (32029.12/180*12.57)
	
	double fReacProb = fSigmaTotalBarnStdr * 2.81865e-3 * 0.6022 / 4.; // bar * gr/cm2 (target areal density.81865 mg/cm2)
	
	//std::cout<<"\n fReacProb: "<<fReacProb<<" fEventCounter: "<<fEventCounter<<std::endl; 
	
	if((fEventCounter-1) %1000 == 0){ 
	
	fReactionEnergyCM = fEbeamCmHist->GetRandom()/1000. * MeV; // MeV	

	fReactionEnergy = fReactionEnergyCM*(fTargetRestMass+fProjectileRestMass)/fTargetRestMass;//MeV

	double rangeBeam = fRangeVsBeamEnergyLeila.Eval(fBeamEnergy / MeV);
	double rangeReaction = fRangeVsBeamEnergyLeila.Eval(fReactionEnergy / MeV);

	fReactionZ = (rangeBeam-rangeReaction);	
	fReactionZ = (fReactionZ * 1000. * TRexSettings::Get()->GetTargetMaterialDensity() / (mg/cm2)*10. - TRexSettings::Get()->GetTargetPhysicalLength()/2.) * mm;
	
	//std::cout<<"\n fReactionZ: "<<fReactionZ<<" fReactionZ/(mg/cm2): "<<fReactionZ/(mg/cm2)<<" fReactionZ*(mg/cm2): "<<fReactionZ*(mg/cm2)<<" reaction: "<<fReacProbA<<" eventno: "<<fEventCounter<<std::endl;
	
	//std::cout<<"\n TRexSettings::Get()->GetTargetMaterialDensity(): "<<TRexSettings::Get()->GetTargetMaterialDensity()<<" TRexSettings::Get()->GetTargetMaterialDensity()/(mg/cm2): "<<TRexSettings::Get()->GetTargetMaterialDensity()/(mg/cm2)<<" TRexSettings::Get()->GetTargetMaterialDensity()*(mg/cm2): "<<TRexSettings::Get()->GetTargetMaterialDensity()*(mg/cm2)<<std::endl;
	
	//std::cout<<"\n TRexSettings::Get()->GetTargetThickness(): "<<TRexSettings::Get()->GetTargetThickness()<<" TRexSettings::Get()->GetTargetThickness()/(mg/cm2): "<<TRexSettings::Get()->GetTargetThickness()/(mg/cm2)<<" TRexSettings::Get()->GetTargetThickness()*(mg/cm2): "<<TRexSettings::Get()->GetTargetThickness()*(mg/cm2)<<std::endl;	
	
    }

    else fReactionEnergyCM = -1.0;*/

}

void TRexBeam::CreateNtupleBranches() {
    // w/ Griffinv10
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    fNtupleID = analysisManager->CreateNtuple("treeGen", "generator output from TRexBeam"); 
    fNtupleColID[0] = analysisManager->CreateNtupleDColumn(fNtupleID, "beamEnergy");
    fNtupleColID[1] = analysisManager->CreateNtupleDColumn(fNtupleID, "beamWidth");
    fNtupleColID[2] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionEnergy");
    fNtupleColID[3] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionEnergyCM");
    fNtupleColID[4] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionX");
    fNtupleColID[5] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionY");
    fNtupleColID[6] = analysisManager->CreateNtupleDColumn(fNtupleID, "reactionZ");
    fNtupleColID[7] = analysisManager->CreateNtupleDColumn(fNtupleID, "thetaCM");
    fNtupleColID[8] = analysisManager->CreateNtupleDColumn(fNtupleID, "ejectileTheta");
    fNtupleColID[9] = analysisManager->CreateNtupleDColumn(fNtupleID, "recoilTheta");
    fNtupleColID[10] = analysisManager->CreateNtupleDColumn(fNtupleID, "ejectilePhi");
    fNtupleColID[11] = analysisManager->CreateNtupleDColumn(fNtupleID, "recoilPhi");
    fNtupleColID[12] = analysisManager->CreateNtupleDColumn(fNtupleID, "ejectileEnergy");
    fNtupleColID[13] = analysisManager->CreateNtupleDColumn(fNtupleID, "recoilEnergy");
    fNtupleColID[14] = analysisManager->CreateNtupleIColumn(fNtupleID, "projectileZ");
    fNtupleColID[15] = analysisManager->CreateNtupleIColumn(fNtupleID, "projectileA");
    fNtupleColID[16] = analysisManager->CreateNtupleIColumn(fNtupleID, "targetZ");
    fNtupleColID[17] = analysisManager->CreateNtupleIColumn(fNtupleID, "targetA");
    fNtupleColID[18] = analysisManager->CreateNtupleIColumn(fNtupleID, "ejectileZ");
    fNtupleColID[19] = analysisManager->CreateNtupleIColumn(fNtupleID, "ejectileA");
    fNtupleColID[20] = analysisManager->CreateNtupleIColumn(fNtupleID, "recoilZ");
    fNtupleColID[21] = analysisManager->CreateNtupleIColumn(fNtupleID, "recoilA");
    fNtupleColID[22] = analysisManager->CreateNtupleDColumn(fNtupleID, "scatteringProbability");
    fNtupleColID[23] = analysisManager->CreateNtupleIColumn(fNtupleID, "reaction");
    fNtupleColID[24] = analysisManager->CreateNtupleDColumn(fNtupleID, "gammaTheta", *fGammaTheta);
    fNtupleColID[25] = analysisManager->CreateNtupleDColumn(fNtupleID, "gammaPhi", *fGammaPhi);
    fNtupleColID[26] = analysisManager->CreateNtupleDColumn(fNtupleID, "gammaEnergy", *fGammaEnergy);
    analysisManager->FinishNtuple(fNtupleID);
    G4cout << "created ntuple treeGen" << G4endl;

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
	fEjectileTheta = fEjectileLab.theta() / radian;
	fEjectilePhi = fEjectileLab.phi() / radian;
	fEjectileEnergy = (fEjectileLab.e() - fEjectileRestMass) / keV;
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
	fRecoilTheta = fRecoilLab.theta() / radian;
	fRecoilPhi = fRecoilLab.phi() / radian;
	fRecoilEnergy = (fRecoilLab.e() - fRecoilRestMass) / keV;
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
		fGammaTheta->push_back((*fGammaLab)[i].theta() / radian);
		fGammaPhi->push_back((*fGammaLab)[i].phi() / radian);
		fGammaEnergy->push_back((*fGammaLab)[i].e() / keV);

		//std::cout << "fGammaEnergy[" << i << "] = " << (*fGammaEnergy)[i] << std::endl;
	}
}

// w/ Griffinv10
void TRexBeam::FillNtuple() { 
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[0], fBeamEnergy);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[1], fBeamWidth);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[2], fReactionEnergy);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[3], fReactionEnergyCM);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[4], fReactionX);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[5], fReactionY);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[6], fReactionZ);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[7], fThetaCM);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[8], fEjectileTheta);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[9], fRecoilTheta);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[10], fEjectilePhi);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[11], fRecoilPhi);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[12], fEjectileEnergy);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[13], fRecoilEnergy);
    analysisManager->FillNtupleIColumn(fNtupleID, fNtupleColID[14], fProjectileZ);
    analysisManager->FillNtupleIColumn(fNtupleID, fNtupleColID[15], fProjectileA);
    analysisManager->FillNtupleIColumn(fNtupleID, fNtupleColID[16], fTargetZ);
    analysisManager->FillNtupleIColumn(fNtupleID, fNtupleColID[17], fTargetA);
    analysisManager->FillNtupleIColumn(fNtupleID, fNtupleColID[18], fEjectileZ);
    analysisManager->FillNtupleIColumn(fNtupleID, fNtupleColID[19], fEjectileA);
    analysisManager->FillNtupleIColumn(fNtupleID, fNtupleColID[20], fRecoilZ);
    analysisManager->FillNtupleIColumn(fNtupleID, fNtupleColID[21], fRecoilA);
    analysisManager->FillNtupleDColumn(fNtupleID, fNtupleColID[22], fScatteringProbability);
    analysisManager->FillNtupleIColumn(fNtupleID, fNtupleColID[23], fReaction);
    analysisManager->AddNtupleRow(fNtupleID);
}

