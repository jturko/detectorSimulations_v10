//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: TRexMessenger.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TRexMessenger.hh"
#include "TRexSettings.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRexMessenger::TRexMessenger(PrimaryGeneratorAction * pgen) :
    fPrimaryGeneratorAction(pgen)
{
    fSetPrimaryGeneratorCmd = new G4UIcmdWithAString("/DetSys/miniball/SetPrimaryGenerator",this);
    fSetPrimaryGeneratorCmd->SetGuidance("set the primary generator type (AngularDistribution, Rutherford, AlphaSource, etc...)");
    fSetPrimaryGeneratorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSimulateEjectilesCmd = new G4UIcmdWithABool("/DetSys/miniball/SimulateEjectiles",this);
    fSimulateEjectilesCmd->SetGuidance("set bool for simulating ejectiles");
    fSimulateEjectilesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSimulateGammasCmd = new G4UIcmdWithABool("/DetSys/miniball/SimulateGammas",this);
    fSimulateGammasCmd->SetGuidance("set bool for simulating gammas");
    fSimulateGammasCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fIncludeEnergyResolutionCmd = new G4UIcmdWithAnInteger("/DetSys/miniball/IncludeEnergyResolution",this);
    fIncludeEnergyResolutionCmd->SetGuidance("set for including an energy resolution (?)");
    fIncludeEnergyResolutionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fIncludeVacuumChamberCmd = new G4UIcmdWithAnInteger("/DetSys/miniball/IncludeVacuumChamber",this);
    fIncludeVacuumChamberCmd->SetGuidance("set for including the vacuum chamber (?)");
    fIncludeVacuumChamberCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSetVacuumChamberTypeCmd = new G4UIcmdWithAString("/DetSys/miniball/SetVacuumChamberType",this);
    fSetVacuumChamberTypeCmd->SetGuidance("set the vacuum chamber type (?)");
    fSetVacuumChamberTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetVacuumChamberGasCmd = new G4UIcmdWithAString("/DetSys/miniball/SetVacuumChamberGas",this);
    fSetVacuumChamberGasCmd->SetGuidance("set the vacuum chamber gas (?)");
    fSetVacuumChamberGasCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
    fSetVacuumChamberGasPressureCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetVacuumChamberGas",this);   
    fSetVacuumChamberGasPressureCmd->SetGuidance("set the vacuum chamber gas pressure");
    fSetVacuumChamberGasPressureCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fSetVacuumChamberGasPressureCmd->SetUnitCategory("Pressure");

    
    fSetTestSourceEnergyCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetTestSourceEnergy",this);
    fSetTestSourceEnergyCmd->SetGuidance("set the test source energy");
    fSetTestSourceEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fSetTestSourceEnergyCmd->SetUnitCategory("Energy");
    
    fSetBeamEnergyCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetBeamEnergy",this);
    fSetBeamEnergyCmd->SetGuidance("set the beam energy");
    fSetBeamEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fSetBeamEnergyCmd->SetUnitCategory("Energy");
    
    fSetBeamWidthCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetBeamWidth",this);
    fSetBeamWidthCmd->SetGuidance("set the beam width");
    fSetBeamWidthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fSetBeamWidthCmd->SetUnitCategory("Length");

    fSetThetaCmMinCmd = new G4UIcmdWithADouble("/DetSys/miniball/SetThetaCmMin",this);
    fSetThetaCmMinCmd->SetGuidance("set the COM minimum theta value (?)");
    fSetThetaCmMinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


    fSetProjectileZCmd = new G4UIcmdWithAnInteger("/DetSys/miniball/SetProjectileZ",this);
    fSetProjectileZCmd->SetGuidance("set the projectile Z");
    fSetProjectileZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetProjectileACmd = new G4UIcmdWithAnInteger("/DetSys/miniball/SetProjectileA",this);
    fSetProjectileACmd->SetGuidance("set the projectile A");
    fSetProjectileACmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetTargetZCmd = new G4UIcmdWithAnInteger("/DetSys/miniball/SetTargetZ",this);
    fSetTargetZCmd->SetGuidance("set the Target Z");
    fSetTargetZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetTargetACmd = new G4UIcmdWithAnInteger("/DetSys/miniball/SetTargetA",this);
    fSetTargetACmd->SetGuidance("set the Target A");
    fSetTargetACmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetEjectileZCmd = new G4UIcmdWithAnInteger("/DetSys/miniball/SetEjectileZ",this);
    fSetEjectileZCmd->SetGuidance("set the ejectile Z");
    fSetEjectileZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetEjectileACmd = new G4UIcmdWithAnInteger("/DetSys/miniball/SetEjectileA",this);
    fSetEjectileACmd->SetGuidance("set the ejectile A");
    fSetEjectileACmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetRecoilZCmd = new G4UIcmdWithAnInteger("/DetSys/miniball/SetRecoilZ",this);
    fSetRecoilZCmd->SetGuidance("set the Recoil Z");
    fSetRecoilZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetRecoilACmd = new G4UIcmdWithAnInteger("/DetSys/miniball/SetRecoilA",this);
    fSetRecoilACmd->SetGuidance("set the Recoil A");
    fSetRecoilACmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    
    fSetProjectileNameCmd = new G4UIcmdWithAString("/DetSys/miniball/SetProjectileName",this);
    fSetProjectileNameCmd->SetGuidance("set the Projectile name");
    fSetProjectileNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetTargetNameCmd = new G4UIcmdWithAString("/DetSys/miniball/SetTargetName",this);
    fSetTargetNameCmd->SetGuidance("set the Target name");
    fSetTargetNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetEjectileNameCmd = new G4UIcmdWithAString("/DetSys/miniball/SetEjectileName",this);
    fSetEjectileNameCmd->SetGuidance("set the Ejectile name");
    fSetEjectileNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetRecoilNameCmd = new G4UIcmdWithAString("/DetSys/miniball/SetRecoilName",this);
    fSetRecoilNameCmd->SetGuidance("set the Recoil name");
    fSetRecoilNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    
    fSetTargetMaterialNameCmd = new G4UIcmdWithAString("/DetSys/miniball/SetTargetMaterialName",this);
    fSetTargetMaterialNameCmd->SetGuidance("set the target material name");
    fSetTargetMaterialNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSetTargetAtomicRatioCmd = new G4UIcmdWithADouble("/DetSys/miniball/SetTargetAtomicRatio",this);
    fSetTargetAtomicRatioCmd->SetGuidance("set the target atomic ratio");
    fSetTargetAtomicRatioCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSetTransferOrCoulexProbabilityCmd = new G4UIcmdWithADouble("/DetSys/miniball/SetTransferOrCoulexProbability",this);
    fSetTransferOrCoulexProbabilityCmd->SetGuidance("set the transfer/coulex reaction probability");
    fSetTransferOrCoulexProbabilityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


    fSetLevelFileCmd = new G4UIcmdWithAString("/DetSys/miniball/SetLevelFile",this);
    fSetLevelFileCmd->SetGuidance("set the location/name of the level file");
    fSetLevelFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetAngularDistributionFileCmd = new G4UIcmdWithAString("/DetSys/miniball/SetAngularDistributionFile",this);
    fSetAngularDistributionFileCmd->SetGuidance("set the location/name of the AngularDistribution file");
    fSetAngularDistributionFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetMassFileCmd = new G4UIcmdWithAString("/DetSys/miniball/SetMassFile",this);
    fSetMassFileCmd->SetGuidance("set the location/name of the Mass file");
    fSetMassFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetCrossSectionFileCmd = new G4UIcmdWithAString("/DetSys/miniball/SetCrossSectionFile",this);
    fSetCrossSectionFileCmd->SetGuidance("set the location/name of the CrossSection file");
    fSetCrossSectionFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetSplineReactionZvsRadiusFileCmd = new G4UIcmdWithAString("/DetSys/miniball/SetSplineReactionZvsRadiusFile",this);
    fSetSplineReactionZvsRadiusFileCmd->SetGuidance("set the location/name of the beam spread file");
    fSetSplineReactionZvsRadiusFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetReactionZDistributionFileCmd = new G4UIcmdWithAString("/DetSys/miniball/SetReactionZDistributionFile",this);
    fSetReactionZDistributionFileCmd->SetGuidance("set the location/name of the reactionZ distribution file");
    fSetReactionZDistributionFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 

    fSetAlphaSourceDiameterCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetAlphaSourceDiameter",this);
    fSetAlphaSourceDiameterCmd->SetGuidance("set the alpha source diameter");
    fSetAlphaSourceDiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fSetAlphaSourceDiameterCmd->SetUnitCategory("Length");
    
    fSetAlphaSourceThicknessCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetAlphaSourceThickness",this);
    fSetAlphaSourceThicknessCmd->SetGuidance("set the alpha source Thickness");
    fSetAlphaSourceThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fSetAlphaSourceThicknessCmd->SetUnitCategory("Length");
    

    fSetTargetDiameterCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetTargetDiameter",this);
    fSetTargetDiameterCmd->SetGuidance("set the target diameter");
    fSetTargetDiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fSetTargetDiameterCmd->SetUnitCategory("Length");
    
    fSetTargetThicknessCmd = new G4UIcmdWithADouble("/DetSys/miniball/SetTargetThickness",this);
    fSetTargetThicknessCmd->SetGuidance("set the target Thickness");
    fSetTargetThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSetGasTargetLengthCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetGasTargetLength",this);
    fSetGasTargetLengthCmd->SetGuidance("set the gas target length");
    fSetGasTargetLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fSetGasTargetLengthCmd->SetUnitCategory("Length");
    
    fSetTargetPressureCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetTargetPressure",this);
    fSetTargetPressureCmd->SetGuidance("set the target pressure");
    fSetTargetPressureCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fSetTargetPressureCmd->SetUnitCategory("Pressure");
    
    
    fSetTargetMaterialDensityCmd = new G4UIcmdWithADouble("/DetSys/miniball/SetTargetMaterialDensity",this);
    fSetTargetMaterialDensityCmd->SetGuidance("set the target material density");
    fSetTargetMaterialDensityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    //fSetTargetMaterialDensityCmd->SetUnitCategory("Volumic Mass"); // this might not be the right unit

    fSetTargetMylarThicknessCmd = new  G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetTargetMylarThickness",this);   
    fSetTargetMylarThicknessCmd->SetGuidance("set the thickness of the mylar foil surrounding the gas target");
    fSetTargetMylarThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSetTargetBeWindowThicknessCmd = new  G4UIcmdWithADoubleAndUnit("/DetSys/miniball/SetTargetBeWindowThickness",this);   
    fSetTargetBeWindowThicknessCmd->SetGuidance("set the thickness of the Be window at the entrance/exit of the gas target");
    fSetTargetBeWindowThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPrintCmd = new G4UIcmdWithoutParameter("/DetSys/miniball/Print",this);
    fPrintCmd->SetGuidance("print values stored in TRexSettings");
    fPrintCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRexMessenger::~TRexMessenger() {
    delete fSetPrimaryGeneratorCmd;
    delete fSimulateEjectilesCmd;
    delete fSimulateGammasCmd;
    delete fIncludeEnergyResolutionCmd;
    delete fIncludeVacuumChamberCmd;
    delete fSetVacuumChamberTypeCmd;
    delete fSetVacuumChamberGasCmd;
    delete fSetVacuumChamberGasPressureCmd;

    delete fSetTestSourceEnergyCmd;
    delete fSetBeamEnergyCmd;
    delete fSetBeamWidthCmd;
    delete fSetThetaCmMinCmd;

    delete fSetProjectileZCmd;
    delete fSetProjectileACmd;
    delete fSetTargetZCmd;
    delete fSetTargetACmd;
    delete fSetEjectileZCmd;
    delete fSetEjectileACmd;
    delete fSetRecoilZCmd;
    delete fSetRecoilACmd;

    delete fSetProjectileNameCmd;
    delete fSetTargetNameCmd;
    delete fSetEjectileNameCmd;
    delete fSetRecoilNameCmd;

    delete fSetTargetMaterialNameCmd;
    delete fSetTargetAtomicRatioCmd;
    delete fSetTransferOrCoulexProbabilityCmd;

    delete fSetLevelFileCmd;
    delete fSetAngularDistributionFileCmd;
    delete fSetMassFileCmd;
    delete fSetCrossSectionFileCmd;
    delete fSetSplineReactionZvsRadiusFileCmd;
    delete fSetReactionZDistributionFileCmd;

    delete fSetAlphaSourceDiameterCmd;
    delete fSetAlphaSourceThicknessCmd;

    delete fSetTargetDiameterCmd;
    delete fSetTargetThicknessCmd;
    delete fSetGasTargetLengthCmd;
    delete fSetTargetPressureCmd;

    delete fSetTargetMaterialDensityCmd;

    delete fSetTargetMylarThicknessCmd;
    delete fSetTargetBeWindowThicknessCmd;

    delete fPrintCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TRexMessenger::SetNewValue(G4UIcommand* command, G4String value) {
    if(command == fSetPrimaryGeneratorCmd)          TRexSettings::Get()->SetPrimaryGenerator(value);
    if(command == fSimulateEjectilesCmd)            TRexSettings::Get()->SimulateEjectiles(fSimulateEjectilesCmd->GetNewBoolValue(value));
    if(command == fSimulateGammasCmd)               TRexSettings::Get()->SimulateGammas(fSimulateGammasCmd->GetNewBoolValue(value));
    if(command == fIncludeEnergyResolutionCmd)      TRexSettings::Get()->IncludeEnergyResolution(fIncludeEnergyResolutionCmd->GetNewIntValue(value));
    if(command == fIncludeVacuumChamberCmd)         TRexSettings::Get()->IncludeVacuumChamber(fIncludeVacuumChamberCmd->GetNewIntValue(value));
    if(command == fSetVacuumChamberTypeCmd)         TRexSettings::Get()->SetVacuumChamberType(value);
    if(command == fSetVacuumChamberGasCmd)          TRexSettings::Get()->SetVacuumChamberGas(value);
    if(command == fSetVacuumChamberGasPressureCmd)  TRexSettings::Get()->SetVacuumChamberGasPressure(fSetVacuumChamberGasPressureCmd->GetNewDoubleValue(value));
    
    if(command == fSetTestSourceEnergyCmd)          TRexSettings::Get()->SetTestSourceEnergy(fSetTestSourceEnergyCmd->GetNewDoubleValue(value));
    if(command == fSetBeamEnergyCmd)                TRexSettings::Get()->SetBeamEnergy(fSetBeamEnergyCmd->GetNewDoubleValue(value));
    if(command == fSetBeamWidthCmd)                 TRexSettings::Get()->SetBeamWidth(fSetBeamWidthCmd->GetNewDoubleValue(value));
    if(command == fSetThetaCmMinCmd)                TRexSettings::Get()->SetThetaCmMin(fSetThetaCmMinCmd->GetNewDoubleValue(value));
    
    if(command == fSetProjectileZCmd)               TRexSettings::Get()->SetProjectileZ(fSetProjectileZCmd->GetNewIntValue(value));
    if(command == fSetProjectileACmd)               TRexSettings::Get()->SetProjectileA(fSetProjectileACmd->GetNewIntValue(value));
    if(command == fSetTargetZCmd)                   TRexSettings::Get()->SetTargetZ(fSetTargetZCmd->GetNewIntValue(value));
    if(command == fSetTargetACmd)                   TRexSettings::Get()->SetTargetA(fSetTargetACmd->GetNewIntValue(value));
    if(command == fSetEjectileZCmd)                 TRexSettings::Get()->SetEjectileZ(fSetEjectileZCmd->GetNewIntValue(value));
    if(command == fSetEjectileACmd)                 TRexSettings::Get()->SetEjectileA(fSetEjectileACmd->GetNewIntValue(value));
    if(command == fSetRecoilZCmd)                   TRexSettings::Get()->SetRecoilZ(fSetRecoilZCmd->GetNewIntValue(value));
    if(command == fSetRecoilACmd)                   TRexSettings::Get()->SetRecoilA(fSetRecoilACmd->GetNewIntValue(value));
    
    if(command == fSetProjectileNameCmd)            TRexSettings::Get()->SetProjectileName(value);
    if(command == fSetTargetNameCmd)                TRexSettings::Get()->SetTargetName(value);
    if(command == fSetEjectileNameCmd)              TRexSettings::Get()->SetEjectileName(value);
    if(command == fSetRecoilNameCmd)                TRexSettings::Get()->SetRecoilName(value);

    if(command == fSetTargetMaterialNameCmd)            TRexSettings::Get()->SetTargetMaterialName(value);
    if(command == fSetTargetAtomicRatioCmd)             TRexSettings::Get()->SetTargetAtomicRatio(fSetTargetAtomicRatioCmd->GetNewDoubleValue(value));
    if(command == fSetTransferOrCoulexProbabilityCmd)   TRexSettings::Get()->SetTransferOrCoulexProbability(fSetTransferOrCoulexProbabilityCmd->GetNewDoubleValue(value));
    
    if(command == fSetLevelFileCmd)                     TRexSettings::Get()->SetLevelFile(value);
    if(command == fSetAngularDistributionFileCmd)       TRexSettings::Get()->SetAngularDistributionFile(value);
    if(command == fSetMassFileCmd)                      TRexSettings::Get()->SetMassFile(value);
    if(command == fSetCrossSectionFileCmd)              TRexSettings::Get()->SetCrossSectionFile(value);
    if(command == fSetSplineReactionZvsRadiusFileCmd)   TRexSettings::Get()->SetSplineReactionZvsRadiusFile(value);
    if(command == fSetReactionZDistributionFileCmd)     TRexSettings::Get()->SetReactionZDistributionFile(value);

    if(command == fSetAlphaSourceDiameterCmd)       TRexSettings::Get()->SetAlphaSourceDiameter(fSetAlphaSourceDiameterCmd->GetNewDoubleValue(value));
    if(command == fSetAlphaSourceThicknessCmd)      TRexSettings::Get()->SetAlphaSourceThickness(fSetAlphaSourceThicknessCmd->GetNewDoubleValue(value));

    if(command == fSetTargetDiameterCmd)            TRexSettings::Get()->SetTargetDiameter(fSetTargetDiameterCmd->GetNewDoubleValue(value));
    if(command == fSetTargetThicknessCmd)           TRexSettings::Get()->SetTargetThickness(fSetTargetThicknessCmd->GetNewDoubleValue(value));
    if(command == fSetGasTargetLengthCmd)           TRexSettings::Get()->SetGasTargetLength(fSetGasTargetLengthCmd->GetNewDoubleValue(value));
    if(command == fSetTargetPressureCmd)            TRexSettings::Get()->SetTargetPressure(fSetTargetPressureCmd->GetNewDoubleValue(value));

    if(command == fSetTargetMaterialDensityCmd)     TRexSettings::Get()->SetTargetMaterialDensity(fSetTargetMaterialDensityCmd->GetNewDoubleValue(value));

    if(command == fSetTargetMylarThicknessCmd)      TRexSettings::Get()->SetTargetMylarThickness(fSetTargetMylarThicknessCmd->GetNewDoubleValue(value));
    if(command == fSetTargetBeWindowThicknessCmd)   TRexSettings::Get()->SetTargetBeWindowThickness(fSetTargetBeWindowThicknessCmd->GetNewDoubleValue(value));

    if(command == fPrintCmd)                        TRexSettings::Get()->Print();

}   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
