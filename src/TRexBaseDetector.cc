/*
 * TRexBaseDetector.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#include "TRexBaseDetector.hh"

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
//=======================================================================
// ctors/dtors
//=======================================================================
TRexBaseDetector::TRexBaseDetector() {
	// TODO Auto-generated constructor stub

}

TRexBaseDetector::~TRexBaseDetector() {
	// TODO Auto-generated destructor stub
}
//=======================================================================
// protected fkt
//=======================================================================
std::vector<ParticleMC>* TRexBaseDetector::CombineDeltaEErestData(TRexBaseDetector * deltaEDet, TRexBaseDetector * restDet) const {
	std::vector<ParticleMC>* particleMCvector = deltaEDet->GetParticleMCvector();
	std::vector<ParticleMC>* particleMCErest = restDet->GetParticleMCvector();

	// nothing to be done if the Erest was not hit
	if(particleMCErest->size() == 0) {
		delete particleMCErest;
		particleMCErest = 0;
		return particleMCvector;
	}

	// fill the data of the Erest detector in the vector of the deltaE detector
	for(unsigned int iDeltaE = 0; iDeltaE < particleMCvector->size(); iDeltaE++) {
		for(unsigned int iErest = iDeltaE; iErest < particleMCErest->size(); iErest++) {
			// check if the detectors have the same id
			if((*particleMCvector)[iDeltaE].GetID() == (*particleMCErest)[iErest].GetID()) {
				(*particleMCvector)[iDeltaE].SetEdet((*particleMCErest)[iErest].GetEdet());
				(*particleMCvector)[iDeltaE].SetStopped((*particleMCErest)[iErest].IsStopped());
				(*particleMCvector)[iDeltaE].SetResKinEnergy((*particleMCErest)[iErest].GetResKinEnergy());
				(*particleMCvector)[iDeltaE].SetA((*particleMCErest)[iErest].GetA());
				(*particleMCvector)[iDeltaE].SetZ((*particleMCErest)[iErest].GetZ());
				(*particleMCvector)[iDeltaE].SetTime((*particleMCErest)[iErest].GetTime());
				(*particleMCvector)[iDeltaE].SetTrackID((*particleMCErest)[iErest].GetTrackID());
			}
		}
	}

	delete particleMCErest;
	particleMCErest = 0;
	return particleMCvector;
}
