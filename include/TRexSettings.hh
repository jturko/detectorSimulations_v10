/*
 * TRexSettings.hh
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#ifndef TREXSETTINGS_HH_
#define TREXSETTINGS_HH_

#include "TEnv.h"

#include "G4SystemOfUnits.hh"

//#include "globals.hh"

#include <string>

class TRexSettings : public TObject 
{
public:
    static TRexSettings* Get();
    virtual ~TRexSettings();
    
    std::string GetPrimaryGenerator() { return fPrimaryGenerator; }
    void SetPrimaryGenerator(std::string pgen) { fPrimaryGenerator = pgen; }
    
    bool SimulateEjectiles() { return fSimulateEjectiles; }
    void SimulateEjectiles(bool val) { fSimulateEjectiles = val; }
    
    bool SimulateGammas() { return fSimulateGammas; }
    void SimulateGammas(bool val) { fSimulateGammas = val; }
    
    int IncludeEnergyResolution() { return fIncludeEnergyResolution; }
    void IncludeEnergyResolution(int val) { fIncludeEnergyResolution = val; }
    
    int IncludeVacuumChamber() { return fIncludeVacuumChamber; }
    void IncludeVacuumChamber(int val) { fIncludeVacuumChamber = val; }
    
    std::string GetVacuumChamberType() { return fVacuumChamberType; }
    void SetVacuumChamberType(std::string type) { fVacuumChamberType = type; }
    
    std::string GetVacuumChamberGas() { return fVacuumChamberGas; }
    void SetVacuumChamberGas(std::string gas) { fVacuumChamberGas = gas; }
    
    double GetVacuumChamberGasPressure() { return fVacuumChamberGasPressure; }
    void SetVacuumChamberGasPressure(double pressure) { fVacuumChamberGasPressure = pressure; }
    
    void Print(Option_t* opt = "") const;
    
    // test source parameters
    double GetTestSourceEnergy() { return fTestSourceEnergy; }
    void SetTestSourceEnergy(double energy) { fTestSourceEnergy = energy; }
    
    // beam parameters
    double GetBeamEnergy() { return fBeamEnergy; }
    void SetBeamEnergy(double val) { fBeamEnergy = val; }
    
    double GetBeamWidth() { return fBeamWidth; }
    void SetBeamWidth(double val) { fBeamWidth = val; }
    
    double GetThetaCmMin() { return fThetaCmMin; }
    void SetThetaCmMin(double val) { fThetaCmMin = val; }
    
    int GetProjectileZ() { return fProjectileZ; }
    int GetProjectileA() { return fProjectileA; }
    int GetTargetZ() { return fTargetZ; }
    int GetTargetA() { return fTargetA; }
    int GetEjectileZ() { return fEjectileZ; }
    int GetEjectileA() { return fEjectileA; }
    int GetRecoilZ() { return fRecoilZ; }
    int GetRecoilA() { return fRecoilA; }
    std::string GetProjectileName() { return fProjectileName; }
    std::string GetTargetName() { return fTargetName; }
    std::string GetEjectileName() { return fEjectileName; }
    std::string GetRecoilName() { return fRecoilName; }
    
    void SetProjectileZ(int val) { fProjectileZ = val; }
    void SetProjectileA(int val) { fProjectileA = val; }
    void SetTargetZ(int val) { fTargetZ = val; }
    void SetTargetA(int val) { fTargetA = val; }
    void SetEjectileZ(int val) { fEjectileZ = val; }
    void SetEjectileA(int val) { fEjectileA = val; }
    void SetRecoilZ(int val) { fRecoilZ = val; }
    void SetRecoilA(int val) { fRecoilA = val; }
    void SetProjectileName(std::string name){ fProjectileName = name; }
    void SetTargetName(std::string name)    { fTargetName = name; }
    void SetEjectileName(std::string name)  { fEjectileName = name; }
    void SetRecoilName(std::string name)    { fRecoilName = name; }
    
    std::string GetTargetMaterialName() { return fTargetMaterialName; }
    void SetTargetMaterialName(std::string name) { fTargetMaterialName = name; }
    
    double GetTargetAtomicRatio() { return fTargetAtomicRatio; }
    void SetTargetAtomicRatio(double ratio) { fTargetAtomicRatio = ratio; }
    
    double GetTransferOrCoulexProbability() { return fTransferOrCoulexProbability; }
    void SetTransferOrCoulexProbability(double prob) { fTransferOrCoulexProbability = prob; }
    
    std::string GetLevelFile() { return fLevelFile; }
    std::string GetAngularDistributionFile() { return fAngularDistributionFile; }
    std::string GetMassFile() { return fMassFile; }
    std::string GetCrossSectionFile() { return fCrossSectionFile; } // Leila
    
    void SetLevelFile(std::string name) { fLevelFile = name; }
    void SetAngularDistributionFile(std::string name) { fAngularDistributionFile = name; }
    void SetMassFile(std::string name) {fMassFile = name; }
    void SetCrossSectionFile(std::string name) { fCrossSectionFile = name; } // Leila
    
    // alpha source  parameters
    double GetAlphaSourceDiameter() { return fAlphaSourceDiameter; }
    void SetAlphaSourceDiameter(double val) { fAlphaSourceDiameter = val; }
    
    double GetAlphaSourceThickness() { return fAlphaSourceThickness; }
    void SetAlphaSourceThickness(double val) { fAlphaSourceThickness = val; }
    
    // target dimensions (or inactive alpha source dimensions)
    double GetTargetDiameter() { return fTargetDiameter; }
    double GetTargetThickness() { return fTargetThickness; }
    double GetTargetThicknessMgPerCm2();
    double GetGasTargetLength() { return fGasTargetLength; }
    double GetTargetPressure() { return fTargetPressure; }
    
    void SetTargetDiameter(double val) { fTargetDiameter =  val; }
    void SetTargetThickness(double val) { fTargetThickness = val*(mg/cm2); }
    void SetGasTargetLength(double val) { fGasTargetLength = val; }
    void SetTargetPressure(double val) { fTargetPressure = val; }
    
    // returns the target length in um
    double GetTargetPhysicalLength();
    
    //std::string GetTargetMaterial() { return fTargetMaterial; }
    double GetTargetMaterialDensity() { return fTargetMaterialDensity; }
    void SetTargetMaterialDensity(double density) { fTargetMaterialDensity = density*(g/cm3); }
    //void SetTargetThickness(double thickness) { fTargetThickness = thickness; }
    
    // target foils
    double GetTargetMylarThickness() { return fTargetMylarThickness; }
    void SetTargetMylarThickness(double thickness) { fTargetMylarThickness = thickness; }

    double GetTargetBeWindowThickness() { return fTargetBeWindowThickness; }
    void SetTargetBeWindowThickness(double thickness) { fTargetBeWindowThickness = thickness; }

    // this constructor has to be public to be able to
    // write the class to file, but it should never be used!
    TRexSettings();

private:
    static TRexSettings* fSettings;
    
    std::string fSettingsFile;
    std::string fPrimaryGenerator;
    bool fSimulateEjectiles;
    bool fSimulateGammas;
    int fIncludeEnergyResolution;
    int fIncludeVacuumChamber;
    std::string fVacuumChamberType;
    std::string fVacuumChamberGas;
    double fVacuumChamberGasPressure;
    
    double fTestSourceEnergy;
    double fBeamEnergy;
    double fBeamWidth;
    double fThetaCmMin;
    
    int fProjectileZ, fProjectileA;
    int fTargetZ, fTargetA;
    int fEjectileZ, fEjectileA;
    int fRecoilZ, fRecoilA;
    std::string fProjectileName, fTargetName, fEjectileName, fRecoilName;
    
    std::string fTargetMaterialName;
    double fTargetAtomicRatio;
    double fTransferOrCoulexProbability;
    
    std::string fLevelFile;
    std::string fAngularDistributionFile;
    std::string fMassFile;
    std::string fCrossSectionFile; // Leila
    
    double fAlphaSourceDiameter;
    double fAlphaSourceThickness;
    
    double fTargetDiameter;
    double fTargetThickness;
    double fGasTargetLength;
    double fTargetPressure;
    //std::string fTargetMaterial;
    double fTargetMaterialDensity;
    
    double fTargetMylarThickness;
    double fTargetBeWindowThickness;

};

#endif /* TREXSETTINGS_HH_ */
