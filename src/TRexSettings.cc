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
}

TRexSettings::~TRexSettings() {
	// TODO Auto-generated destructor stub
}

void TRexSettings::ReadSettingsFile(std::string settingsFile) {
	fSettingsFile = settingsFile; 
	TEnv sett(fSettingsFile.c_str());

	fPrimaryGenerator = sett.GetValue("PrimaryGenerator", "beam");

	fSimulateEjectiles = sett.GetValue("SimulateEjectiles", false);
	fSimulateGammas = sett.GetValue("SimulateGammas", false);
	fIncludeEnergyResolution = sett.GetValue("IncludeEnergyResolution", 0);

	//vacuum chamber
	fIncludeVacuumChamber = sett.GetValue("IncludeVacuumChamber", 1);
	fVacuumChamberType = sett.GetValue("VacuumChamberType", "cylinder");
	fVacuumChamberGas = sett.GetValue("VacuumChamberGas", "helium");
	//fVacuumChamberGasPressure = sett.GetValue("VacuumChamberGasPressure", 1e-6) /1000. * bar; // original
	//fVacuumChamberGasPressure = sett.GetValue("VacuumChamberGasPressure", 0.) /1000. * bar; // no gas pressure in the chamber
	fVacuumChamberGasPressure = sett.GetValue("VacuumChamberGasPressure", sett.GetValue("TargetPressure", 1000.0) /1000.) * bar;

	fTestSourceEnergy = sett.GetValue("TestSourceEnergy", 5000.0) * keV;

	//beam properties
	fBeamEnergy = sett.GetValue("BeamEnergy", 100.0) * MeV;
	fBeamWidth = sett.GetValue("BeamWidth", 3.0) * mm; //original
	//fBeamWidth = RandFlat::shoot(sett.GetValue("BeamWidth", 3.0) * (-0.5) * mm, sett.GetValue("BeamWidth", 3.0) * 0.5 * mm);
	//fThetaCmMin = sett.GetValue("ThetaCmMin", 2.0) * degree; // commented out by Leila 27.07.2017
	fThetaCmMin = sett.GetValue("ThetaCmMin", 2.0) * degree;  // added by Leila 27.07.2017
	
	//G4cout<<"Leilaaaaaaaaaaaaaaaaaaaaaaaaaaa RexSettings: fThetaCmMin "<<fThetaCmMin<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555"<<G4endl;
	//G4cout<<"Leilaaaaaaaaaaaaaaaaaaaaaaaaaaa RexSettings: vacuumChamberGasPressure "<<fVacuumChamberGasPressure/bar <<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555"<<G4endl;

	//reaction
	fProjectileZ = sett.GetValue("ProjectileZ", 30);
	fProjectileA = sett.GetValue("ProjectileA", 72);
	fTargetZ = sett.GetValue("TargetZ", 1);
	fTargetA = sett.GetValue("TargetA", 3);
	fEjectileZ = sett.GetValue("EjectileZ", 30);
	fEjectileA = sett.GetValue("EjectileA", 74);
	fRecoilZ = sett.GetValue("RecoilZ", 1);
	fRecoilA = sett.GetValue("RecoilA", 1);
	fProjectileName = sett.GetValue("ProjectileName", "72Zn");
	fTargetName = sett.GetValue("TargetName", "3H");
	fEjectileName = sett.GetValue("EjectileName", "74Zn");
	fRecoilName = sett.GetValue("RecoilName", "1H");

	std::string simDir = "./";
	fLevelFile = sett.GetValue("LevelFile", "dummy.dat");
	if(fLevelFile.find("LevelFiles/") == 0) {
		fLevelFile.insert(0, "/");
		fLevelFile.insert(0, simDir);
	}
	fAngularDistributionFile = sett.GetValue("AngularDistributionFile", "dummy.dat");
	if(fAngularDistributionFile.find("AngularDistributionFiles/") == 0) {
		fAngularDistributionFile.insert(0, "/");
		fAngularDistributionFile.insert(0, simDir);
	}
	fCrossSectionFile = sett.GetValue("CrossSectionFile", "dummy.dat");
	if(fCrossSectionFile.find("AngularDistributionFiles/") == 0) {
		fCrossSectionFile.insert(0, "/");
		fCrossSectionFile.insert(0, simDir);
	}
	fMassFile = sett.GetValue("MassFile", "MassFile.dat");
	if(fMassFile.find("MassFile.dat") == 0) {
		fMassFile.insert(0, "/");
		fMassFile.insert(0, simDir);
	}

	fAlphaSourceDiameter = sett.GetValue("AlphaSourceDiameter", 3.) * mm; // original 3 mm //last value: 2
	fAlphaSourceThickness = sett.GetValue("AlphaSourceThickness",3.) * mm;// original 3 mm //last value 40 

	fTargetDiameter = sett.GetValue("TargetDiameter", 3) * mm;
	fTargetThickness = sett.GetValue("TargetThickness", 0.5) * mg/cm2;
	fGasTargetLength = sett.GetValue("GasTargetLength", 0.0) * cm;
	fTargetPressure = sett.GetValue("TargetPressure", 0.0) /1000. * bar;
	//fTargetMaterial = sett.GetValue("TargetMaterial", "vacuum");
	fTargetMaterialDensity = sett.GetValue("TargetMaterialDensity", 4.507) * g/cm3;

	fTargetMaterialName = sett.GetValue("TargetMaterialName", "dummy");
	//std::cout<<"from '"<<fTargetMaterialName<<"' to ";
	//std::transform(fTargetMaterialName.begin(), fTargetMaterialName.end(), fTargetMaterialName.begin(), ::tolower);
	//std::cout<<"'"<<fTargetMaterialName<<"'"<<std::endl;
	fTargetAtomicRatio = sett.GetValue("TargetAtomicRatio", 1.5);
	fTransferOrCoulexProbability = sett.GetValue("TransferOrCoulexProbability", 1.0);

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
