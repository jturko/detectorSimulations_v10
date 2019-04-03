/*
 * TRexRutherford.hh
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 * 
 * Modified 2017/06/15 dhymers
 * Added isDefined to keep track of first-run instantiations
 */


#ifndef TREXRUTHERFORD_HH_
#define TREXRUTHERFORD_HH_

#include "TRexBeam.hh"
#include "G4ParticleGun.hh"

#include "Randomize.hh"
#include "TMath.h"


class TRexRutherford : public TRexBeam {
	public:
		TRexRutherford();
		virtual ~TRexRutherford();

		void GeneratePrimaries(G4Event *anEvent);

	private:
		void ShootEjectileAndRecoilDirections();
		void CalculateScatteringProbability();
		void DoKinematicCalculation();

		G4double fNorm;
		G4double fEcm;
		G4double fTcm;
		G4double fTcm3;
		G4double fTcm4;
		G4double fPcm3;
		G4double fPcm4;
		G4double fBeta;
		G4double fGamma;
		G4double fBetaCm3;
		G4double fBetaCm4;
		
		G4bool isDefined = false;
};

#endif /* TREXRUTHERFORD_HH_ */
