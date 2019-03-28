/*
 * TRexBaseGenerator.hh
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#ifndef TREXBASEGENERATOR_HH_
#define TREXBASEGENERATOR_HH_

#include "TTree.h"

#include "G4Event.hh"

class TRexBaseGenerator {
	public:
		TRexBaseGenerator();
		virtual ~TRexBaseGenerator();

		//virtual void GeneratePrimaries(G4Event *anEvent) {};
		virtual void GeneratePrimaries(G4Event*) {}
		//virtual void CreateTreeBranches(TTree &tree) {};
		//virtual void FillTree(TTree &tree);
		virtual void CreateTreeBranches() {}
		virtual void FillTree();

		void SetTree(TTree* tree);

	protected:
		TTree* fTree;
};

#endif /* TREXBASEGENERATOR_HH_ */
