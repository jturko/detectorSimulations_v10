/*
 * TRexAlphaSource.hh
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

//Simulates a quadruple alpha source

#ifndef TREXALPHASOURCE_HH_
#define TREXALPHASOURCE_HH_

#include "TRexBaseGenerator.hh"

//#include "math.h"

class G4ParticleGun;

class TRexAlphaSource : public TRexBaseGenerator {
	public:
		TRexAlphaSource();
		virtual ~TRexAlphaSource();

		void GeneratePrimaries(G4Event *anEvent);
		//void CreateTreeBranches(TTree &tree);
		//void FillTree(TTree &tree);
		void CreateNtupleBranches(TTree * tree);
        void FillNtuple();

	private:
		// shoot fReactionX, fReactionY, fReactionZ at the surface of the alpha source
		void ShootReactionPosition();

		// shoot fThetaCM and fPhi isotropically in a half sphere
		void CreateIsotropicDistribution();

		G4ParticleGun* fParticleGun;

		G4double fReactionEnergy;
		G4double fReactionX;
		G4double fReactionY;
		G4double fReactionZ;
		G4double fThetaCM, fPhi;
};

#endif /* TREXALPHASOURCE_HH_ */
