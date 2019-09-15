/*
 * TRexBaseDetector.hh
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#ifndef TREXBASEDETECTOR_HH_
#define TREXBASEDETECTOR_HH_

#include "ParticleMC.hh"

#include "G4UnitsTable.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"

class G4VSolid;
//class G4LogicalVolume;
class G4VPhysicalVolume;

#include <vector>
#include <string>

class TRexBaseDetector {
	public:
		TRexBaseDetector();
		TRexBaseDetector(std::string name, std::string direction, int id);
		virtual ~TRexBaseDetector();

		//virtual void Construct(G4LogicalVolume* experimentalHall_log, G4SDManager *SDMan) {};
		virtual void Build(G4LogicalVolume*, G4SDManager*) {}
		//virtual void CreateTreeBranches() = 0;
		//virtual ParticleMC* GetParticleMC() = 0;
		virtual std::vector<ParticleMC>* GetParticleMCvector() = 0;

		std::string GetName() {return fName; }


	protected:
		std::string fName;
		std::string fBaseName;
		std::string fDirection;
		int fId;

		G4VSolid* fSolid;
		G4LogicalVolume* fLogicalVolume;
		G4VPhysicalVolume* fPhysicalVolume;

		//! combines the data from a DeltaE and a rest detector
		//! @param deltaEDet the deltaE detector
		//! @param restDet the rest detector
		//! @return the combined particle MC vector, the caller takes ownership of the vector
		virtual std::vector<ParticleMC>* CombineDeltaEErestData(TRexBaseDetector * deltaEDet, TRexBaseDetector * restDet) const;
};

#endif /* TREXBASEDETECTOR_HH_ */
