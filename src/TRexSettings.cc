/*
 * TRexSettings.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#include "TRexSettings.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
//#include "Randomize.hh"

#include <algorithm>

TRexSettings* TRexSettings::fSettings = NULL;

TRexSettings* TRexSettings::Get() {
	if(fSettings == NULL) {
		fSettings = new TRexSettings;
	}

	return fSettings;
}

TRexSettings::TRexSettings() 
{
    fPrimaryGenerator = "beam";
    
    fSimulateEjectiles = false;
    fSimulateGammas = false;
    fIncludeEnergyResolution = 0;

    fIncludeVacuumChamber = 1;
    fVacuumChamberType = "cylinder";
    fVacuumChamberGas = "helium";
    fVacuumChamberGasPressure = 1.0;

    fTestSourceEnergy = 5000.*keV;
    
    fBeamEnergy = 100.*MeV;
    fBeamWidth = 3.*mm;
    fThetaCmMin = 2.*degree;
	
    fProjectileZ = 30;
	fProjectileA = 72;
	fTargetZ = 1;
	fTargetA = 3;
	fEjectileZ = 30;
	fEjectileA = 74;
	fRecoilZ = 1;
	fRecoilA = 1;
	fProjectileName = "72Zn";
	fTargetName = "3H";
	fEjectileName = "74Zn";
	fRecoilName = "1H";

    fLevelFile = "";
    fAngularDistributionFile = "";
    fCrossSectionFile = "";
    fMassFile = "";
    fReactionZvsRadiusFile = "";
    fReactionZDistributionFile = "";	

    fAlphaSourceDiameter = 3. * mm; // original 3 mm //last value: 2
	fAlphaSourceThickness = 3. * mm;// original 3 mm //last value 40 

	fTargetDiameter = 3.0 * mm;
	fTargetThickness = 0.5 * mg/cm2;
	fGasTargetLength = 0.0 * cm;
	fTargetPressure = 0./1000.0 * bar;
	fTargetMaterialDensity = 4.507 * g/cm3;

	fTargetMaterialName = "dummy";
	fTargetAtomicRatio = 1.5;
	fTransferOrCoulexProbability = 1.0;

    fTargetMylarThickness = 2.*um;
    fTargetBeWindowThickness = 8.*um;

    fReactionZvsRadiusFileBool = false;
    fReactionZDistributionFileBool = false;
}

TRexSettings::~TRexSettings() {
	// TODO Auto-generated destructor stub
}

void TRexSettings::Print(Option_t* opt) const {
	std::cout<<"TRexSettings: "<<opt<<std::endl
		<<"fPrimaryGenerator = "<<fPrimaryGenerator<<std::endl
		<<"fSimulateEjectiles = "<<fSimulateEjectiles<<std::endl
		<<"fSimulateGammas = "<<fSimulateGammas<<std::endl
		<<"fIncludeEnergyResolution = "<<fIncludeEnergyResolution<<std::endl
		<<"fIncludeVacuumChamber = "<<fIncludeVacuumChamber<<std::endl
		<<"fVacuumChamberType = "<<fVacuumChamberType<<std::endl
		<<"fVacuumChamberGas = "<<fVacuumChamberGas<<std::endl
		<<"fVacuumChamberGasPressure = "<<fVacuumChamberGasPressure*1000./bar<<" mbar"<<std::endl
		<<"fTestSourceEnergy = "<<fTestSourceEnergy/MeV<<" MeV"<<std::endl
		<<"---------- beam properties"<<std::endl
		<<"fBeamEnergy = "<<fBeamEnergy/MeV<<" MeV"<<std::endl
		<<"fBeamWidth = "<<fBeamWidth/mm<<" mm"<<std::endl
		<<"fThetaCmMin = "<<fThetaCmMin/degree<<" degree"<<std::endl
		<<"---------- reaction"<<std::endl
		<<"fProjectileZ = "<<fProjectileZ<<std::endl
		<<"fProjectileA = "<<fProjectileA<<std::endl
		<<"fTargetZ = "<<fTargetZ<<std::endl
		<<"fTargetA = "<<fTargetA<<std::endl
		<<"fEjectileZ = "<<fEjectileZ<<std::endl
		<<"fEjectileA = "<<fEjectileA<<std::endl
		<<"fRecoilZ = "<<fRecoilZ<<std::endl
		<<"fRecoilA = "<<fRecoilA<<std::endl
		<<"fProjectileName = "<<fProjectileName<<std::endl
		<<"fTargetName = "<<fTargetName<<std::endl
		<<"fEjectileName = "<<fEjectileName<<std::endl
		<<"fRecoilName = "<<fRecoilName<<std::endl
		<<"---------- other files"<<std::endl
		<<"fLevelFile = "<<fLevelFile<<std::endl
		<<"fAngularDistributionFile = "<<fAngularDistributionFile<<std::endl
		<<"fCrossSectionFile = "<<fCrossSectionFile<<std::endl
		<<"fMassFile = "<<fMassFile<<std::endl
		<<"fReactionZvsRadiusFile = "<<fReactionZvsRadiusFile<<std::endl
		<<"fReactionZDistributionFile = "<<fReactionZDistributionFile<<std::endl
		<<"---------- file booleans"<<std::endl
		<<"fReactionZvsRadiusFileBool = "<<fReactionZvsRadiusFileBool<<std::endl
		<<"fReactionZDistributionFileBool = "<<fReactionZDistributionFileBool<<std::endl
		<<"---------- alpha source"<<std::endl
		<<"fAlphaSourceDiameter = "<<fAlphaSourceDiameter/mm<<" mm"<<std::endl
		<<"fAlphaSourceThickness = "<<fAlphaSourceThickness/mm<<" mm"<<std::endl
		<<"---------- target"<<std::endl
		<<"fTargetDiameter = "<<fTargetDiameter/mm<<" mm"<<std::endl
		<<"fTargetThickness = "<<fTargetThickness/(mg/cm2)<<" mg/cm2"<<std::endl
		<<"fGasTargetLength = "<<fGasTargetLength/mm<<" mm"<<std::endl
		<<"fTargetPressure = "<<fTargetPressure*1000./bar<<" mbar"<<std::endl
		<<"fTargetMaterialDensity = "<<fTargetMaterialDensity/(g/cm3)<<" g/cm3"<<std::endl
		<<"fTargetMaterialName = "<<fTargetMaterialName<<std::endl
		<<"fTargetAtomicRatio = "<<fTargetAtomicRatio<<std::endl
        <<"fTargetMylarThickness = "<<fTargetMylarThickness/um<<" um"<<std::endl
        <<"fTargetBeWindowThickness = "<<fTargetBeWindowThickness/um<<" um"<<std::endl
		<<"fTransferOrCoulexProbability = "<<fTransferOrCoulexProbability<<std::endl
		<<std::endl;
}

G4double TRexSettings::GetTargetPhysicalLength(){
	if(GetGasTargetLength() > 0.) {
		return GetGasTargetLength();
	}
	G4double thickness = GetTargetThickness() / GetTargetMaterialDensity() ;

	return (thickness);
}

double TRexSettings::GetTargetThicknessMgPerCm2() { 
	return fTargetThickness/(mg/cm2); 
}
