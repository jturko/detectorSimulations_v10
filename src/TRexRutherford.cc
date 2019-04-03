/*
 * TRexRutherford.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 *
 * Modified 2017/06/15 trockman
 * Call DefineNuclei after physics list is instantiated
 * 
 * Modified 2017/06/15 dhymers
 * Moved additional things into first-run call to GeneratePrimaries to
 * fix some dependency issues.
 * 
 * fThetaCmMin = sett.GetValue("ThetaCmMin", 2.0) * CLHEP::degree; from TRexSettings.cc
 * 
 * 
 */

#include "TRexRutherford.hh"
#include "TRexSettings.hh"

TRexRutherford::TRexRutherford() {
	// normalization constant
	fNorm = 1./  (1./(sin(fThetaCM_min/CLHEP::rad * 0.5) * sin(fThetaCM_min/CLHEP::rad * 0.5)) -1);

	// calculate reaction energy in the middle of the target
	//fReactionZ = 0.;
	//CalculateReactionEnergyInTheTarget();

	// calculate the probability for a Rutherford scattering
	//CalculateScatteringProbability();

	// kinematics calculation before the reaction
	//DoKinematicCalculation();
}

TRexRutherford::~TRexRutherford() {
	// TODO Auto-generated destructor stub
}


void TRexRutherford::GeneratePrimaries(G4Event *anEvent) {
	if (isDefined == false){
		//define nuclei after physics list is instantiated
		DefineNuclei();
		
		fTargetMaterial = GetTargetMaterial();
		std::cout << "TargetMaterialName for energy loss calculation in the target = " << fTargetMaterial->Name() << std::endl;
		fKinematics = new Kinematic(&fProjectile, fTargetMaterial, TRexSettings::Get()->GetTargetThickness()/(CLHEP::mg/CLHEP::cm2));
		
		fEnergyVsTargetDepth = *(fKinematics->EnergyVsThickness(fBeamEnergy / CLHEP::MeV, TRexSettings::Get()->GetTargetThickness() / 1000 / (CLHEP::mg/CLHEP::cm2)));
		
		// calculate reaction energy in the middle of the target
		fReactionZ = 0.;
		CalculateReactionEnergyInTheTarget();
		
		// calculate the probability for a Rutherford scattering
		CalculateScatteringProbability();
		
		// kinematics calculation before the reaction
		DoKinematicCalculation();
		
		isDefined = true;
	}

	// shoot the emission point
	ShootReactionPosition();

	// calculate reaction energy in the target
	CalculateReactionEnergyInTheTarget();

	// shoot energy and direction
	ShootEjectileAndRecoilDirections();

	// set the gun properties and create the primary vertex
	SetEjectileGun(anEvent);
	SetRecoilGun(anEvent);

	FillNtuple();
}


void TRexRutherford::CalculateScatteringProbability() {
	// Rutherford factor
	G4double RF = fProjectileZ * fTargetZ * CLHEP::eplus * CLHEP::eplus / (16. * M_PI * CLHEP::epsilon0);
	RF *= RF;
	//G4cout<<"RF = "<<RF/(MeV * MeV *millibarn)<<G4endl;

	// total energy in the CM frame
	//G4double E_CM = projectileCM.e() - projectileRestMass +
	//  targetCM.e() - targetRestMass; DAS IST FALSCH!!!!
	G4double E_CM = fReactionEnergy * fTargetRestMass / (fProjectileRestMass + fTargetRestMass);

	G4double sigmaTot = -4. * M_PI * RF / (E_CM * E_CM) * (1. - 1./ (sin(fThetaCM_min *0.5) * sin(fThetaCM_min * 0.5)));
	//G4cout<<"sigmaTot = "<<sigmaTot/millibarn<<G4endl;

	G4double arealDensity = TRexSettings::Get()->GetTargetThickness() * CLHEP::Avogadro / (fTargetA* CLHEP::g/CLHEP::mole);
	//G4cout<<"areal density = "<<arealDensity*cm2<<G4endl;

	fScatteringProbability = arealDensity * sigmaTot;
	//G4cout<<"scatteringProbability = "<<fScatteringProb<<G4endl;
}


void TRexRutherford::DoKinematicCalculation() {
	// kinematics (RelativistischeKinematik.pdf)
	fEcm = sqrt(pow(fReactionEnergy + fProjectileRestMass + fTargetRestMass,2) -
			pow(fReactionEnergy,2) - 2. * fReactionEnergy * fProjectileRestMass);

	//std::cout<<"fReactionEnergy = "<<fReactionEnergy / MeV<<" , fBeamEnergy = "<<fBeamEnergy / MeV<<" , fEcm = "<<fEcm / MeV<<std::endl;

	fTcm = fEcm - fProjectileRestMass - fTargetRestMass;

	fTcm3 = 0.5 * fTcm * (fTcm + 2. * fTargetRestMass) / fEcm;
	fTcm4 = 0.5 * fTcm * (fTcm + 2. * fProjectileRestMass) / fEcm;

	fPcm3 = sqrt(fTcm3 * fTcm3 + 2. * fTcm3 * fProjectileRestMass) / CLHEP::c_light;
	fPcm4 = sqrt(fTcm4 * fTcm4 + 2. * fTcm4 * fTargetRestMass) / CLHEP::c_light;

	fBeta = sqrt(pow(fReactionEnergy, 2) + 2. * fProjectileRestMass * fReactionEnergy) / (fReactionEnergy + fProjectileRestMass + fTargetRestMass);
	fGamma = 1. / sqrt(1. - fBeta * fBeta);

	fBetaCm3 = fPcm3 * CLHEP::c_light / (fProjectileRestMass + fTcm3);
	fBetaCm4 = fPcm4 * CLHEP::c_light / (fTargetRestMass + fTcm4);
}


void TRexRutherford::ShootEjectileAndRecoilDirections() {
	// particle momentum direction
	fEjectilePhi = CLHEP::RandFlat::shoot(-M_PI, M_PI) * CLHEP::rad;
	fRecoilPhi = -fEjectilePhi;

	// shoot theta according Rutherford differential cross section
	G4double rand = CLHEP::RandFlat::shoot(0., 1.);
	/* not correct as the sin(thetaCM) from the spherical coordinates are missing
	  do {
	    thetaCM = CLHEP::RandFlat::shoot(thetaCM_min/CLHEP::rad,M_PI);
	    rand = CLHEP::RandFlat::shoot(0., RF / (E_CM * E_CM) / pow(sin(thetaCM_min/2./CLHEP::rad),4));
	  } while(rand > RF / (E_CM * E_CM) / pow(sin(thetaCM/2./CLHEP::rad),4));
	 */
	fThetaCM = 2.* asin(1. / sqrt(1. / sin(fThetaCM_min*0.5) / sin(fThetaCM_min*0.5) - rand / fNorm));

	G4ThreeVector ejectileMomentumVectorCM;
	ejectileMomentumVectorCM.setRThetaPhi(fPcm3, fThetaCM, fEjectilePhi);

	G4LorentzVector ejectileCM;
	ejectileCM.setE(fProjectileRestMass + fTcm3);
	ejectileCM.setVect(ejectileMomentumVectorCM);

	G4ThreeVector recoilMomentumVectorCM;
	recoilMomentumVectorCM.setRThetaPhi(fPcm4, TMath::Pi() - fThetaCM, TMath::Pi() + fEjectilePhi);

	G4LorentzVector recoilCM;
	recoilCM.setE(fTargetRestMass + fTcm4);
	recoilCM.setVect(recoilMomentumVectorCM);

	// transformation from the CM frame to the lab frame
	double T3 = (fGamma - 1.) * fProjectileRestMass + fGamma * fTcm3 + fGamma *fBeta * fPcm3 * CLHEP::c_light * cos(fThetaCM);
	double T4 = (fGamma - 1.) * fTargetRestMass + fGamma * fTcm4 + fGamma * fBeta * fPcm4 * CLHEP::c_light * cos(TMath::Pi() - fThetaCM);

	double p3 = sqrt(T3 * T3 + 2. * T3 * fProjectileRestMass) / CLHEP::c_light;
	double p4 = sqrt(T4 * T4 + 2. * T4 * fTargetRestMass) / CLHEP::c_light;

	double theta3 = fabs(atan(sin(fThetaCM) / (fGamma * (cos(fThetaCM) + fBeta / fBetaCm3))));
	double theta4 = fabs(atan(sin(TMath::Pi() - fThetaCM) / (fGamma * (cos(TMath::Pi() - fThetaCM) + fBeta / fBetaCm4))));

	// construct Lorentz vectors in the lab frame
	G4ThreeVector ejectileMomentumVectorLab;
	ejectileMomentumVectorLab.setRThetaPhi(p3, theta3, fEjectilePhi);

	//std::cout<<"theta3 = "<<theta3 * 180. /  TMath::Pi()<<std::endl;
	//std::cout<<"theta4 = "<<theta4 * 180. /  TMath::Pi()<<std::endl;

	fEjectileLab.setE(fProjectileRestMass + T3);
	fEjectileLab.setVect(ejectileMomentumVectorLab);

	G4ThreeVector recoilMomentumVectorLab;
	recoilMomentumVectorLab.setRThetaPhi(p4, theta4, TMath::Pi() + fEjectilePhi);

	fRecoilLab.setE(fTargetRestMass + T4);
	fRecoilLab.setVect(recoilMomentumVectorLab);
}
