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
// $Id: TRexMessenger.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TREXMESSENGER_HH
#define TREXMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

class TRexSettings;

class PrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TRexMessenger: public G4UImessenger
{
public:
    TRexMessenger(PrimaryGeneratorAction* pgen);
    virtual ~TRexMessenger();

public:
    void SetNewValue(G4UIcommand*, G4String);

private:

    PrimaryGeneratorAction*     fPrimaryGeneratorAction;

    G4UIcmdWithAString*         fSetPrimaryGeneratorCmd;
    G4UIcmdWithABool*           fSimulateEjectilesCmd;
    G4UIcmdWithABool*           fSimulateGammasCmd;
    G4UIcmdWithAnInteger*       fIncludeEnergyResolutionCmd;
    G4UIcmdWithAnInteger*       fIncludeVacuumChamberCmd;
    G4UIcmdWithAString*         fSetVacuumChamberTypeCmd;
    G4UIcmdWithAString*         fSetVacuumChamberGasCmd;
    G4UIcmdWithADoubleAndUnit*  fSetVacuumChamberGasPressureCmd;

    G4UIcmdWithADoubleAndUnit*  fSetTestSourceEnergyCmd;
    G4UIcmdWithADoubleAndUnit*  fSetBeamEnergyCmd;
    G4UIcmdWithADoubleAndUnit*  fSetBeamWidthCmd;
    G4UIcmdWithADouble*  fSetThetaCmMinCmd;

    G4UIcmdWithAnInteger* fSetProjectileZCmd;
    G4UIcmdWithAnInteger* fSetProjectileACmd;
    G4UIcmdWithAnInteger* fSetTargetZCmd;
    G4UIcmdWithAnInteger* fSetTargetACmd;
    G4UIcmdWithAnInteger* fSetEjectileZCmd;
    G4UIcmdWithAnInteger* fSetEjectileACmd;
    G4UIcmdWithAnInteger* fSetRecoilZCmd;
    G4UIcmdWithAnInteger* fSetRecoilACmd;

    G4UIcmdWithAString* fSetProjectileNameCmd;
    G4UIcmdWithAString* fSetTargetNameCmd;
    G4UIcmdWithAString* fSetEjectileNameCmd;
    G4UIcmdWithAString* fSetRecoilNameCmd;

    G4UIcmdWithAString*     fSetTargetMaterialNameCmd;
    G4UIcmdWithADouble*     fSetTargetAtomicRatioCmd;
    G4UIcmdWithADouble*     fSetTransferOrCoulexProbabilityCmd;
    
    G4UIcmdWithAString*     fSetLevelFileCmd;
    G4UIcmdWithAString*     fSetAngularDistributionFileCmd;
    G4UIcmdWithAString*     fSetMassFileCmd;
    G4UIcmdWithAString*     fSetCrossSectionFileCmd;

    G4UIcmdWithADoubleAndUnit*  fSetAlphaSourceDiameterCmd;
    G4UIcmdWithADoubleAndUnit*  fSetAlphaSourceThicknessCmd;

    G4UIcmdWithADoubleAndUnit*  fSetTargetDiameterCmd;
    G4UIcmdWithADouble*         fSetTargetThicknessCmd;
    G4UIcmdWithADoubleAndUnit*  fSetGasTargetLengthCmd;
    G4UIcmdWithADoubleAndUnit*  fSetTargetPressureCmd;

    G4UIcmdWithADouble*         fSetTargetMaterialDensityCmd;
 
    G4UIcmdWithoutParameter*    fPrintCmd;
   
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

