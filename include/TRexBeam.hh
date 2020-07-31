/*
 * TRexBeam.hh
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

// Simulates a beam

#ifndef TREXBEAMSOURCE_HH_
#define TREXBEAMSOURCE_HH_

#include "TRexBaseGenerator.hh"
#include "Element.hh"
#include "Isotopes.hh"
#include "Kinematic.hh"

#include "G4ParticleTable.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Gamma.hh"
#include "G4LorentzVector.hh"

#include "Randomize.hh"

#include "math.h"

#include "TSpline.h"
//#include "TRandom.h"
//#include "TRandom1.h"
#include "TH1F.h"

#include "TFile.h"

class G4ParticleGun;

class TRexBeam : public TRexBaseGenerator 
{
public:
    TRexBeam();
    virtual ~TRexBeam();
    
    virtual void GeneratePrimaries(G4Event*) {}
    //void GeneratePrimaries(G4Event *anEvent);
    //void CreateTreeBranches(TTree &tree);
    //void FillTree(TTree &tree);
    void FillNtuple();
    void CreateNtupleBranches(TTree * tree);
	
    void SaveExtras(TFile * file);
    virtual double GetReactionRatio() {}
    
protected:
    // shoot fReactionX, fReactionY, fReactionZ of the source
    void ShootReactionPosition();
    
    // build a spline that defines the reaction x-y radius as a function of z, the depth in the target
    void BuildReactionZvsRadiusSpline();

    // fill 
    void FillReactionZDistributionHisto();
    void FillReactionZDistributionGraph();

    // define nuclei
    void DefineNuclei();
    
    // to calculate energy loss of the beam in the target
    Material* GetTargetMaterial();
    void FillCrossSectionGraph(); // Leila
    void CalculateReactionEnergyInTheTarget();
    
    // define particle type
    G4ParticleDefinition* ParticleDefinition(int Z, int N, double eex);
    
    // set particle guns
    void SetEjectileGun(G4Event *anEvent);
    void SetRecoilGun(G4Event *anEvent);
    void SetGammaGun(G4Event *anEvent);
    
    G4ParticleGun* fParticleGunEjectile;
    G4ParticleGun* fParticleGunRecoil;
    G4ParticleGun* fParticleGunGamma;
    
    Isotopes* fIsotopeTable;
    
    G4int fProjectileZ;
    G4int fProjectileA;
    G4int fTargetZ;
    G4int fTargetA;
    G4int fEjectileZ;
    G4int fEjectileA;
    G4int fRecoilZ;
    G4int fRecoilA;
    
    G4double fProjectileRestMass;
    G4double fTargetRestMass;
    G4double fEjectileRestMass;
    G4double fRecoilRestMass;
    
    Material* fTargetMaterial;
    Kinematic* fKinematics;
    Element fProjectile;
    Element fTarget;
    Element fEjectile;
    Element fRecoil;
    TSpline3 fEnergyVsTargetDepth;
    TSpline3 fRangeVsBeamEnergyLeila;
    
    G4double fBeamEnergy;
    G4double fBeamWidth;
    G4double fReactionEnergy;
    G4double fReactionEnergyCM;
    G4double fReactionX;
    G4double fReactionY;
    G4double fReactionZ;
    G4double fThetaCM;
    G4double fEjectileTheta;
    G4double fEjectilePhi;
    G4double fEjectileEnergy;
    G4double fRecoilTheta;
    G4double fRecoilPhi;
    G4double fRecoilEnergy;
    std::vector<G4double>* fGammaTheta;
    std::vector<G4double>* fGammaPhi;
    std::vector<G4double>* fGammaEnergy;
    
    G4double fScatteringProbability;
    G4double fThetaCM_min;
    size_t fReaction;
    
    G4double fExcitationEnergy;
    
    G4LorentzVector fEjectileLab;
    G4LorentzVector fRecoilLab;
    G4LorentzVector fEjectileLabAfterTarget;
    std::vector<G4LorentzVector> *fGammaLab;
    
    int fNbOfBeamEnergyInCm; // Leila
    std::vector<double> fEbeamCm; // Leila
    std::vector<double> fsigmaForEbeamCm; // Leila
    std::vector<TGraph> fGraphCrossSection;// Leila 
    G4double fSigmaVsEbeamCmMax;//Leila
    TGraph* fGrp;
    //TH1F* fEbeamCmHist;
    int fNumberOfPointsGraph;
    double* fYaxs;
    //TRandom fRndReaction;
    int fEventCounter;

    TSpline3 * fReactionZvsRadiusSpline;

    TGraph * fReactionZDistributionGraph;
    TH1F * fReactionZDistributionHisto;
    TSpline3 * fReactionZDistributionSpline;
};

#endif /* TREXBEAMSOURCE_HH_ */
