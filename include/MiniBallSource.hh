/************************************************************************
 * \file MiniBallSource.hh
 *
 * \brief mother class for referrence sources
 *
 * \author hans.boie@mpi-hd.mpg.de
 *
 ************************************************************************/

#ifndef MiniBallSource_H
#define MiniBallSource_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include <vector>

#define UNKNOWN -1
#define E0 0
#define E1 1
#define E2 2
#define M1 3
#define M1_E2 10
#define INT_CONV 99

#define MB_SOURCE_DL 0

/************************************************************************
 * \class MiniBallSourceLevel
 *
 * \brief energy level of source
 *
 ************************************************************************/
class MiniBallSourceLevel
{
public:
  MiniBallSourceLevel() {;};
  MiniBallSourceLevel(G4double en) { energy = en;};
  MiniBallSourceLevel(G4double p_en, G4double p_l, G4int p_p)
                                   { energy = p_en; l = p_l; p = p_p; };
  ~MiniBallSourceLevel() {;};

public:

  G4double energy;
  G4double l;
  G4int p;

  std::vector<G4int> decay;
  std::vector<G4int> decay_type;
  std::vector<G4double> delta;
  std::vector<G4double> prob;

  G4long pop_stat;
  G4long feed_stat;
  std::vector<G4long> decay_stat;
  std::vector<G4long> conv_stat;

public:

  void AddDecay(G4int level_id, G4double prob,
		G4int decay_type=UNKNOWN, G4double delta=0.);
};

/************************************************************************
 * \class MiniBallSource
 *
 * \brief possible decay branch for angular correlation
 *
 ************************************************************************/
class MiniBallSourceDecayBranch
{
public:
  MiniBallSourceDecayBranch(G4double m_i, G4double m_f,
			    G4double m, G4double m1, G4double m2,
			    G4double prob);
  ~MiniBallSourceDecayBranch() {;};

public:
  G4double m_i, m_f, m, m1, m2, prob;
};

/************************************************************************
 * \class MiniBallSourceLevel
 *
 * \brief mother class for referrence sources
 *
 ************************************************************************/
class MiniBallSource : public G4VUserPrimaryGeneratorAction
{

public:
  MiniBallSource(G4double x_pos=0., G4double y_pos=0., G4double z_pos=0.);
  ~MiniBallSource();

  enum Parity
    {
      PLUS = +1,
      MINUS = -1
    };
  
  enum SourceType
    {
      AlphaSource,
      BetaMinusSource,
      BetaPlusSource,
      GammaSource
    };

public:
  /** \brief generation of primary particles
   * 
   *  This method will be use if the class is used as PrimaryGeneratorAction
   *  of a simulation. This method also illustrates how you can include
   *  the simulation of the source into your own PrimaryGeneratorAction class.
   */
  void GeneratePrimaries(G4Event*);

  /** \brief include simulation of angular distribution
   * 
   *  For simple decays the simulation of the angular distribution
   *  can be included, e.g. this is used in the class MiniBallCobalt60Source.
   *  Take care that only a few possible cascades are implemented.
   *  Therefore the simulation of angular distribution is truned off 
   *  by default.
   */
  void SimulateAngularDistribution() { angular_distribution=true; };

  /** \brief generation of a gamma cascade
   * 
   *  A gamma cascade is generated and stored. The number of gammas and
   *  the energies and momentum directions are return by the methods
   *  GetNbOfGammas, GetGammaDirection, GetGammaEnergy.
   *
   *  \see GetNbOfGammas
   *  \see GetGammaDirection
   *  \see GetGammaEnergy
   */
  void GenerateDecay();

  /** \brief generation of a gamma cascade from specified energy level
   * 
   *  A gamma cascade from the specified level is generated and stored. 
   *  The number of gammas and
   *  the energies and momentum directions are return by the methods
   *  GetNbOfGammas, GetGammaDirection, GetGammaEnergy.
   *
   *  \see GetNbOfGammas
   *  \see GetGammaDirection
   *  \see GetGammaEnergy
   */
  void GenerateDecay(G4int);

  /** \brief returns the number of gammas in the cascade
   * 
   *  \return number of gammas in the cascade
   */
  G4int GetNbOfGammas();

  /** \brief returns the direction of a given gamma
   * 
   *  If the angular distribution is not simulated, then direction
   *  of the gammas will be distributed isotropicly.
   *
   *  \param number of gamma in the cascade 
   *         (ranges from 0 the GetNbOfGammas()-1)
   *  \return direction of a given gamma
   */
  G4ThreeVector GetGammaDirection(G4int gamma_nb);

  /** \brief returns the direction of a given gamma
   * 
   *  If the angular distribution is not simulated, then direction
   *  of the gammas will be distributed isotropicly.
   *
   *  \param number of gamma in the cascade 
   *         (ranges from 0 the GetNbOfGammas()-1)
   *  \return direction of a given gamma
   */
  G4double GetGammaEnergy(G4int gamma_nb);

  G4double RandTheta(G4int l, G4int m);
  G4double GetF(G4int l, G4int m, G4double theta);
  G4double GetFmax(G4int l, G4int m);
  G4double GetClebschGordan(G4double j1, G4double m1,
			    G4double j2, G4double m2,
			    G4double J, G4double M);

  void SetInternalConversionCoeffE1(G4double c, G4double en);
  void SetInternalConversionCoeffE2(G4double c, G4double en);
  G4double GetInternalConversionCoeffE1(G4double en);
  G4double GetInternalConversionCoeffE2(G4double en);
  G4double GetInternalConversionCoeff(G4int type, G4double en);

  void InternalConversionCorrectProbabilities();

  void ClearStatistics();
  void PrintStatistics();

protected:
  SourceType type;

  G4double source_x, source_y, source_z;

  G4int DecayLevel(G4int level_id);

  G4double GetAngularDistribution(G4int level_id, G4int new_level_id,
				  G4int decay_type);
  G4int GetDecayType(G4int decay_type, G4double delta);
  void GetDecayValues(G4int decay_type,
		      G4double* Ldiff, G4bool* Pdiff);

  void CheckSettings();

public:
  G4ParticleGun* particleGun;

  std::vector<MiniBallSourceLevel*> level_structure;
  MiniBallSourceLevel* feeding;

  std::vector<G4double> gamma_energy;
  std::vector<G4ThreeVector> gamma_direction;
  std::vector<G4RotationMatrix> gamma_transform;
  std::vector<G4int> decay_level_history;
  std::vector<G4int> decay_type_history;

  std::vector<G4double> conv_coeff_E1;
  std::vector<G4double> conv_en_E1;
  std::vector<G4double> conv_coeff_E2;
  std::vector<G4double> conv_en_E2;

  G4bool cleared_statistics;
  G4bool checked_settings;

  G4bool angular_distribution;

  G4long events;

  G4bool ang_err_msg;
};

#endif
