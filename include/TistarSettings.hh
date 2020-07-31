/*
 * TistarSettings.hh
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#ifndef TISTARSETTINGS_HH_
#define TISTARSETTINGS_HH_

#include "TEnv.h"
#include "TVector3.h"

#include "G4SystemOfUnits.hh"

//#include "globals.hh"

#include <iostream>
#include <string>

class TistarSettings : public TObject 
{
public:
    static TistarSettings* Get();
    virtual ~TistarSettings();
    
    static void Set(TistarSettings*);

    std::string GetPrimaryGenerator() { return fPrimaryGenerator; }
    std::string GetSecondPrimaryGenerator() { return fSecondPrimaryGenerator; }

    void SetPrimaryGenerator(std::string pgen) { fPrimaryGenerator = pgen; }
    void SetSecondPrimaryGenerator(std::string spgen) { fSecondPrimaryGenerator = spgen; }

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
   
    std::string GetVacuumChamberMaterialName() { return fVacuumChamberMaterialName; }
    void SetVacuumChamberMaterialName(std::string name) { fVacuumChamberMaterialName = name; }

    double GetVacuumChamberGasPressure() { return fVacuumChamberGasPressure; }
    void SetVacuumChamberGasPressure(double pressure) { fVacuumChamberGasPressure = pressure; }

    double GetVacuumChamberGasTemperature() { return fVacuumChamberGasTemperature; }
    void SetVacuumChamberGasTemperature(double temperature) { fVacuumChamberGasTemperature = temperature; }

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
    double GetTargetTemperature() { return fTargetTemperature; }
    
    void SetTargetDiameter(double val) { fTargetDiameter =  val; }
    void SetTargetThickness(double val) { fTargetThickness = val * (CLHEP::mg/CLHEP::cm2); }
    void SetGasTargetLength(double val) { fGasTargetLength = val; }
    void SetTargetPressure(double val) { fTargetPressure = val; }
    void SetTargetTemperature(double val) { fTargetTemperature = val; }    

    // returns the target length in um
    double GetTargetPhysicalLength();
    
    //std::string GetTargetMaterial() { return fTargetMaterial; }
    double GetTargetMaterialDensity() { return fTargetMaterialDensity; }
    void SetTargetMaterialDensity(double density) { fTargetMaterialDensity = density * (CLHEP::g/CLHEP::cm3); }
    //void SetTargetThickness(double thickness) { fTargetThickness = thickness; }
    
    // target foils
    double GetTargetMylarThickness() { return fTargetMylarThickness; }
    void SetTargetMylarThickness(double thickness) { fTargetMylarThickness = thickness; }

    double GetTargetBeWindowThickness() { return fTargetBeWindowThickness; }
    void SetTargetBeWindowThickness(double thickness) { fTargetBeWindowThickness = thickness; }

    // for beam spreading (beam radius) as a function of z
    std::string GetReactionZvsRadiusFile() { return fReactionZvsRadiusFile; }
    void SetReactionZvsRadiusFile(std::string name) { fReactionZvsRadiusFileBool = true; fReactionZvsRadiusFile = name; }
    bool GetReactionZvsRadiusFileBool() { return fReactionZvsRadiusFileBool; }
    
    std::string GetReactionZDistributionFile() { return fReactionZDistributionFile; }
    void SetReactionZDistributionFile(std::string name) { fReactionZDistributionFileBool = true; fReactionZDistributionFile = name; }
    bool GetReactionZDistributionFileBool() { return fReactionZDistributionFileBool; }

    // detectors
    std::vector<double> GetLayerPositionZ(int layerN) { return fLayerPositionZ[layerN]; }
    std::vector<double> GetLayerDistToBeam(int layerN) { return fLayerDistToBeam[layerN]; }
    double GetLayerDimensionsX(int layerN) { return fLayerDimensionsX[layerN]; }
    double GetLayerDimensionsY(int layerN) { return fLayerDimensionsY[layerN]; }
    
    void SetLayerPositionZ(int layerN, std::vector<double> pos) { fLayerPositionZ.resize(layerN+1); fLayerPositionZ[layerN] = pos; }
    void SetLayerDistToBeam(int layerN, std::vector<double> dist) { fLayerDistToBeam.resize(layerN+1); fLayerDistToBeam[layerN] = dist; }
    void SetLayerDimensionsX(int layerN, double dimX) { fLayerDimensionsX.resize(layerN+1); fLayerDimensionsX[layerN] = dimX; }
    void SetLayerDimensionsY(int layerN, double dimY) { fLayerDimensionsY.resize(layerN+1); fLayerDimensionsY[layerN] = dimY; }

    std::vector<std::vector<TVector3>> GetLayerPositionVector() { return fLayerPositionVector; }
    std::vector<std::vector<TVector3>> GetLayerDimensionVector() { return fLayerDimensionVector; }

    void SetLayerPositionVector(int layerN, int stripN, TVector3 vector) { 
        fLayerPositionVector.resize(layerN+1); 
        fLayerPositionVector[layerN].resize(stripN+1); 
        fLayerPositionVector[layerN][stripN] = vector;
    }
    void SetLayerDimensionVector(int layerN, int stripN, TVector3 vector) { 
        fLayerDimensionVector.resize(layerN+1); 
        fLayerDimensionVector[layerN].resize(stripN+1); 
        fLayerDimensionVector[layerN][stripN] = vector;
    }

    // this constructor has to be public to be able to
    // write the class to file, but it should never be used!
    TistarSettings();
    
    void SaveMe(bool val) { fSaveMe = val; }
    bool SaveMe() { return fSaveMe; }

    void Nevents(Int_t n) { fNevents = n; }
    Int_t Nevents() { return fNevents; }

private:
    static TistarSettings* fSettings;
    
    std::string fSettingsFile;
    std::string fPrimaryGenerator;
    std::string fSecondPrimaryGenerator;
    bool fSimulateEjectiles;
    bool fSimulateGammas;
    int fIncludeEnergyResolution;
    int fIncludeVacuumChamber;
    std::string fVacuumChamberType;
    std::string fVacuumChamberGas;
    std::string fVacuumChamberMaterialName;
    double fVacuumChamberGasPressure;
    double fVacuumChamberGasTemperature;

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
    std::string fReactionZvsRadiusFile;    
    std::string fReactionZDistributionFile;

    double fAlphaSourceDiameter;
    double fAlphaSourceThickness;
    
    double fTargetDiameter;
    double fTargetThickness;
    double fGasTargetLength;
    double fTargetPressure;
    double fTargetMaterialDensity;
    double fTargetTemperature;    

    double fTargetMylarThickness;
    double fTargetBeWindowThickness;

    bool fReactionZvsRadiusFileBool;

    bool fReactionZDistributionFileBool;
 
    // detectors
    std::vector<std::vector<double>> fLayerPositionZ;
    std::vector<std::vector<double>> fLayerDistToBeam;
    std::vector<double> fLayerDimensionsX;
    std::vector<double> fLayerDimensionsY;

    std::vector<std::vector<TVector3>> fLayerDimensionVector;
    std::vector<std::vector<TVector3>> fLayerPositionVector;

    // extra
    bool fSaveMe;
    Int_t fNevents;

    ClassDef(TistarSettings, 1);

};

#endif /* TISTARSETTINGS_HH_ */
