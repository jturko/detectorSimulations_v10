/***************************************************************************
 *
 *   Materials used for the miniball detectors
 *
 ***************************************************************************/

#include "MiniBallMaterial.hh"

#include "G4Material.hh"
#include "G4Element.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

MiniBallMaterial* MiniBallMaterial::fMiniBallMaterial = NULL;

MiniBallMaterial* MiniBallMaterial::Get() {
	if(fMiniBallMaterial == NULL) {
		fMiniBallMaterial = new MiniBallMaterial;
	}

	return fMiniBallMaterial;
}

MiniBallMaterial::MiniBallMaterial()
{
	// germanium
	Ge = new G4Material("HPGermanium", 32., 72.61*g/mole, 5.32*g/cm3);

	// aluminum
	Al = new G4Material("Aluminum", 13., 26.981539*g/mole, 2.70*g/cm3);

	//Vacuum 10^-7 bar and room temperature
	N = new G4Element("Nitrogen", "N", 7., 14.01*g/mole);

	O = new G4Element("Oxygen", "O", 8., 16.00*g/mole);

	G4double density = 1.290*mg/cm3;
	Air = new G4Material("Air", density, 2);
	Air->AddElement(N, 0.7);
	Air->AddElement(O, 0.3);

	G4double vac_density = 1.225e-7*mg/cm3; 
	G4double pressure    = 1.e-7*bar;
	G4double temperature = CLHEP::STP_Temperature + 20*kelvin;
	Vacuum = new G4Material("Vacuum_7", vac_density, 1,
			kStateGas, temperature, pressure);
	Vacuum->AddMaterial(Air, 1.); 
}

MiniBallMaterial::~MiniBallMaterial() 
{
	delete Ge;
	delete Al;
	delete N;
	delete O;
	delete Air;
	delete Vacuum;
}

