// Includes Physical Constants and System of Units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//#include "G4FieldManager.hh"
//#include "G4UniformMagField.hh"
//#include "MagneticField.hh"
//#include "G4TransportationManager.hh"
//#include "Field.hh"
//#include "GlobalField.hh"

#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"

#include "DetectionSystem8pi.hh"
#include "DetectionSystemGriffin.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

#include "TistarSettings.hh"

void DetectorConstruction::DefineSuppressedParameters()
{
	fGriffinFwdBackPosition = 11.0*cm;
	fDetectorRadialDistance = 11.0*cm ;
}

void DetectorConstruction::DefineMaterials()
{ 
	// use G4-NIST materials data base
	//
	G4NistManager* man = G4NistManager::Instance();
	man->FindOrBuildMaterial("G4_Galactic");
	man->FindOrBuildMaterial("G4_Pb");
	man->FindOrBuildMaterial("G4_lAr");
	man->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    man->FindOrBuildMaterial("G4_Xe");

	man->FindOrBuildMaterial("G4_Al");
	man->FindOrBuildMaterial("G4_POLYETHYLENE");
	man->FindOrBuildMaterial("G4_RUBBER_NEOPRENE");
	man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	man->FindOrBuildMaterial("G4_BGO");
	man->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	man->FindOrBuildMaterial("G4_Ge");
	man->FindOrBuildMaterial("G4_Cu");
	man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	man->FindOrBuildMaterial("G4_AIR");
	G4Material* SODIUM_IODIDE   = man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
	G4Material* tl              = man->FindOrBuildMaterial("G4_Tl");

	G4double a, z, density, temperature, pressure;
	G4String name, symbol;
	G4int    nElements, nComponents, nAtoms;
	G4double fractionMass;

	std::vector<G4Element*>  myElements;    // save pointers here to avoid
	std::vector<G4Material*> myMaterials;   // warnings of unused variables

	// Elements
	G4Element* elH = new G4Element(name="H", symbol="H", z=1., a = 1.00*g/mole);
	myElements.push_back(elH);

	G4Element* elC = new G4Element(name="C", symbol="C", z=6., a = 12.01*g/mole);
	myElements.push_back(elC);

	G4Element* elN  = new G4Element(name="N", symbol="N",  z=7.,  a= 14.00674*g/mole);
	myElements.push_back(elN);

	G4Element* elO  = new G4Element(name="O",   symbol="O",  z=8.,  a= 15.9994 *g/mole);
	myElements.push_back(elO);

	G4Element* elAr  = new G4Element(name="Ar",   symbol="Ar",  z=18.,  a= 39.948 *g/mole);
	myElements.push_back(elAr);

	G4Element* elLa = new G4Element(name="La", symbol="La", z=57., 138.9055*g/mole);
	myElements.push_back(elLa);

	G4Element* elBr = new G4Element(name="Br", symbol="Br", z=35., 79.904*g/mole);
	myElements.push_back(elBr);

	G4Element* elGe = new G4Element(name="Ge", symbol="Ge", z=32., 72.64*g/mole);
	myElements.push_back(elGe);

	G4Element* elI = new G4Element(name="I", symbol="I", z=53., 126.90*g/mole);
	myElements.push_back(elI);

	G4Element* elCs = new G4Element(name="Cs", symbol="Cs", z=55., 132.91*g/mole);
	myElements.push_back(elCs);

	G4Element* elTa = new G4Element(name="Ta", symbol="Ta", z=73., 180.95*g/mole);
	myElements.push_back(elTa);

	G4Element* elW = new G4Element(name="W", symbol="W", z=74., 183.84*g/mole);
	myElements.push_back(elW);

	G4Element* elBi = new G4Element(name="Bi", symbol="Bi", z=83., 208.98*g/mole);
	myElements.push_back(elBi);

	G4Element* elCe = new G4Element(name="Ce", symbol="Ce", z=58., 140.116*g/mole);
	myElements.push_back(elCe);

	G4Element* elNi = new G4Element(name="Ni", symbol="Ni", z=28., 58.69*g/mole);
	myElements.push_back(elNi);

	G4Element* elCu = new G4Element(name="Cu", symbol="Cu", z=29., 63.546*g/mole);
	myElements.push_back(elCu);

	G4Element* elNd = new G4Element(name="Nd", symbol="Nd", z=60., 144.242*g/mole);
	myElements.push_back(elNd);

	G4Element* elFe = new G4Element(name="Fe", symbol="Fe", z=26., 55.845*g/mole);
	myElements.push_back(elFe);

	G4Element* elB = new G4Element(name="B", symbol="B", z=5., 10.811*g/mole);
	myElements.push_back(elB);

	G4Element* elNa = new G4Element(name="Na", symbol="Na", z=11., 22.99*g/mole);
	myElements.push_back(elNa);

	G4Element* elLi = new G4Element(name="Li", symbol="Li", z=3. , a=  6.94   *g/mole);
	myElements.push_back(elLi);

	G4Element* elF  = new G4Element(name="F" , symbol="F" , z=9. , a= 18.9984 *g/mole);
	myElements.push_back(elF);

	G4Element* elSi = new G4Element(name="Si", symbol="Si", z=14., a= 28.0855 *g/mole);
	myElements.push_back(elSi);

	G4Element* elCr = new G4Element(name="Cr", symbol="Cr", z=24., a= 51.9961 *g/mole);
	myElements.push_back(elCr);

	G4Element* elPb = new G4Element(name="Pb", symbol="Pb", z=82., a=207.2    *g/mole);
	myElements.push_back(elPb);

	G4Element* elPd = new G4Element(name="Palladium", symbol="Pd", z=46., a=106.42    *g/mole);
	myElements.push_back(elPd);

	// Materials
	//  density     = universe_mean_density; //from PhysicalConstants.h
	//  pressure    = 1.e-19*pascal;
	//  temperature = 0.1*kelvin;
	//  G4Material* Galactic = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density, kStateGas,temperature,pressure);
	//  myMaterials.push_back(Galactic);

	G4Material* air = new G4Material(name="Air", density=1.29*mg/cm3, nElements=2);
	air->AddElement(elN, .7);
	air->AddElement(elO, .3);
	myMaterials.push_back(air);

	//  G4Material* TechVacuum = new G4Material( "TechVacuum", density=1.e-5*g/cm3, 1, kStateGas, temperature=STP_Temperature, pressure=2.e-2*bar );
	//  TechVacuum->AddMaterial(air, 1.);
	//  myMaterials.push_back(TechVacuum);

	//  G4Material* Vacuum = new G4Material(name="Vacuum", z=1., a= 1.01*g/mole, density= universe_mean_density, kStateGas, temperature = 0.1*kelvin, pressure=1.0e-19*pascal);
	//  myMaterials.push_back(Vacuum); // "Galatic" Vacuum

	//universe_mean_density = 1.e-25*g/cm3;

	density     = 1.e-24*g/cm3;
	G4Material* vacuumDensityYoctogramPerCm3 = new G4Material(name="vacuumDensityYoctogramPerCm3", z=1., a=1.01*g/mole, density);
	myMaterials.push_back(vacuumDensityYoctogramPerCm3);

	density     = 1.e-9*g/cm3;
	G4Material* vacuumDensityNanogramPerCm3 = new G4Material(name="vacuumDensityNanogramPerCm3", z=1., a=1.01*g/mole, density);
	myMaterials.push_back(vacuumDensityNanogramPerCm3);

	density     = 10.e-9*g/cm3;
	G4Material* vacuumDensity10nanogramPerCm3 = new G4Material(name="vacuumDensity10nanogramPerCm3", z=1., a=1.01*g/mole, density);
	myMaterials.push_back(vacuumDensity10nanogramPerCm3);

	density     = 100.e-9*g/cm3;
	G4Material* vacuumDensity100nanogramPerCm3 = new G4Material(name="vacuumDensity100nanogramPerCm3", z=1., a=1.01*g/mole, density);
	myMaterials.push_back(vacuumDensity100nanogramPerCm3);

	density     = 1.e-6*g/cm3;
	G4Material* vacuumDensityMicrogramPerCm3 = new G4Material(name="vacuumDensityMicrogramPerCm3", z=1., a=1.01*g/mole, density);
	myMaterials.push_back(vacuumDensityMicrogramPerCm3);

	density     = 1.e-3*g/cm3;
	G4Material* vacuumDensityMilligramPerCm3 = new G4Material(name="vacuumDensityMilligramPerCm3", z=1., a=1.01*g/mole, density);
	myMaterials.push_back(vacuumDensityMilligramPerCm3);

	density     = 1.e-2*g/cm3;
	G4Material* vacuumDensityCentigramPerCm3 = new G4Material(name="vacuumDensityCentigramPerCm3", z=1., a=1.01*g/mole, density);
	myMaterials.push_back(vacuumDensityCentigramPerCm3);

	density     = 1.e-1*g/cm3;
	G4Material* vacuumDensityDecigramPerCm3 = new G4Material(name="vacuumDensityDecigramPerCm3", z=1., a=1.01*g/mole, density);
	myMaterials.push_back(vacuumDensityDecigramPerCm3);

	density     = 1.0*g/cm3;
	G4Material* vacuumDensityGramPerCm3 = new G4Material(name="vacuumDensityGramPerCm3", z=1., a=1.01*g/mole, density);
	myMaterials.push_back(vacuumDensityGramPerCm3);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// From http://geant4.cern.ch/support/source/geant4/examples/advanced/purgingMagnet/src/PurgMagDetectorConstruction.cc
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Laboratory vacuum: Dry air (average composition)
	density = 1.7836*mg/cm3 ;       // STP
	G4Material* argon = new G4Material(name="argon", density, nComponents=1);
	argon->AddElement(elAr, 1);

	density = 1.25053*mg/cm3 ;       // STP
	G4Material* nitrogen = new G4Material(name="N2", density, nComponents=1);
	nitrogen->AddElement(elN, 2);

	density = 1.4289*mg/cm3 ;       // STP
	G4Material* oxygen = new G4Material(name="O2", density, nComponents=1);
	oxygen->AddElement(elO, 2);

	// LaboratoryVacuum
	density  = 1.2928*mg/cm3 ;       // STP
	density *= 1.0e-8 ;              // pumped vacuum
	temperature = STP_Temperature;
	pressure = 1.0e-8*STP_Pressure;
	G4Material* vacuum = new G4Material(name="Vacuum",	density, nComponents=3,
			kStateGas, temperature, pressure);
	vacuum->AddMaterial(nitrogen, fractionMass = 0.7557 ) ;
	vacuum->AddMaterial(oxygen,   fractionMass = 0.2315 ) ;
	vacuum->AddMaterial(argon,    fractionMass = 0.0128 ) ;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	G4Material* water = new G4Material(name="water", density=1000*kg/m3, nElements=2);
	water->AddElement(elH, nAtoms=2);
	water->AddElement(elO, nAtoms=1);
	myMaterials.push_back(water);

	G4Material* al = new G4Material(name="Aluminum", z=13., a= 26.98154*g/mole, density= 2.70  *g/cm3);
	myMaterials.push_back(al);

	G4Material* w = new G4Material(name="Tungsten", z=74., a= 183.84*g/mole, density= 19.25*g/cm3);
	myMaterials.push_back(w);

	G4Material* si = new G4Material(name="Silicon", z=14., a= 28.0855*g/mole, density= 2.330  *g/cm3);
	myMaterials.push_back(si);

	G4Material* ti = new G4Material(name="Titanium", z=22., 47.867*g/mole, 4.54*g/cm3); //
	myMaterials.push_back(ti);

	G4Material* sn = new G4Material(name="Tin", z=50., 118.71*g/mole, 6.99*g/cm3);//
	myMaterials.push_back(sn);

	G4Material* au = new G4Material(name="Gold", z=79., a= 196.9666*g/mole, density= 19.30  *g/cm3);
	myMaterials.push_back(au);

	G4Material* pb = new G4Material(name="Lead", z=82., a= 207.19*g/mole, density= 11.35  *g/cm3);
	myMaterials.push_back(pb);

	G4Material* ge = new G4Material(name="Germanium", z=32., a= 72.64*g/mole, density=5.323 *g/cm3);
	myMaterials.push_back(ge);

	G4Material* laBr3 = new G4Material(name="LanthanumBromide", density=5.08*g/cm3, nElements=2);
	laBr3->AddElement(elLa, nAtoms=1);
	laBr3->AddElement(elBr, nAtoms=3);
	myMaterials.push_back(laBr3);

	G4Material* laBrCe = new G4Material(name="CeriumDopedLanthanumBromide", density=5.08*g/cm3, nElements=2);
	laBrCe->AddMaterial(laBr3, 95*perCent);
	laBrCe->AddElement(elCe, 5*perCent);
	myMaterials.push_back(laBrCe);

	G4Material* hevimetal = new G4Material("Hevimetal", density=19.0*g/cm3, nElements=3);
	hevimetal->AddElement(elTa, 80*perCent);
	hevimetal->AddElement(elNi, 13*perCent);
	hevimetal->AddElement(elCu, 7*perCent);
	myMaterials.push_back(hevimetal);

	G4Material* lN2 = new G4Material("LiquidN2", density=804*kg/m3, nComponents=1);
	lN2->AddElement(elN, nAtoms=2);
	myMaterials.push_back(lN2);

	G4Material* delrin = new G4Material("Delrin", density = 1.4*g/cm3, nComponents=3);
	delrin->AddElement(elC, nAtoms=1);
	delrin->AddElement(elO, nAtoms=1);
	delrin->AddElement(elH, nAtoms=2);
	myMaterials.push_back(delrin);

	G4Material* be = new G4Material("Beryllium", z=4., a=9.012*g/mole, density=1848*kg/m3);
	myMaterials.push_back(be);

	G4Material* cu = new G4Material("Copper", z=29., a = 63.546*g/mole, density = 8960*kg/m3);
	myMaterials.push_back(cu);

	G4Material* zn = new G4Material("Zinc", 30., 65.409*g/mole, 7.14*g/cm3);
	myMaterials.push_back(zn);

	G4Material* brass= new G4Material("Brass", 8.5*g/cm3, 2);
	brass->AddMaterial(cu, 70*perCent);
	brass->AddMaterial(zn, 30*perCent);
	myMaterials.push_back(brass);

	G4Material* bgo = new G4Material("BGO", density = 7.3*g/cm3, nComponents=3);
	bgo->AddElement(elBi, nAtoms=4);
	bgo->AddElement(elGe, nAtoms=3);
	bgo->AddElement(elO, nAtoms=12);
	myMaterials.push_back(bgo);

	G4Material* csI = new G4Material("CesiumIodide", density = 4.51*g/cm3, nComponents=2);
	csI->AddElement(elCs, nAtoms=1);
	csI->AddElement(elI, nAtoms=1);
	myMaterials.push_back(csI);

	G4Material* bc404 = new G4Material("BC404", density = 1.032*g/cm3, nComponents=2);
	bc404->AddElement(elH, 52.4*perCent);
	bc404->AddElement(elC, 47.6*perCent);
	myMaterials.push_back(bc404);

	G4Material* mylar = new G4Material("Mylar", density = 1.397*g/cm3, nComponents=3);
	mylar->AddElement(elC, nAtoms=10);
	mylar->AddElement(elH, nAtoms=8);
	mylar->AddElement(elO, nAtoms=4);
	myMaterials.push_back(mylar);

	G4Material* ndFeB = new G4Material("NdFeB", density = 7.45*g/cm3, nComponents=3);  //From Wikipedia
	ndFeB->AddElement(elNd, nAtoms=2);
	ndFeB->AddElement(elFe, nAtoms=14);
	ndFeB->AddElement(elB, nAtoms=1);
	myMaterials.push_back(ndFeB);

	G4Material* teflon = new G4Material("Teflon", density = 2.2*g/cm3, nComponents=2);
	teflon->AddElement(elC, nAtoms=2);
	teflon->AddElement(elF, nAtoms=4);
	myMaterials.push_back(teflon);

	G4Material* siLi = new G4Material("Si(Li)", density = 2.330*g/cm3, nComponents=2);
	siLi->AddElement(elSi, 99.999*perCent);
	siLi->AddElement(elLi,  0.001*perCent);
	myMaterials.push_back(siLi);

	G4Material* plate = new G4Material("Plate", density = 55.84*g/cm3, nComponents=3);
	plate->AddElement(elFe, 90.*perCent);
	plate->AddElement(elCr, 9. *perCent);
	plate->AddElement(elPb, 1. *perCent);
	myMaterials.push_back(plate);

	G4Material* naI = new G4Material("NaI", density = 3.67*g/cm3, nComponents=2);
	naI->AddElement(elNa, nAtoms=1);
	naI->AddElement(elI,  nAtoms=1);
	naI->GetIonisation()->SetMeanExcitationEnergy(452*eV);
	myMaterials.push_back(naI);

	G4Material* g4SodiumIodidePlusTl = new G4Material("G4SODIUMIODIDEPLUSTl", density=3.667*g/cm3, nComponents = 2);
	g4SodiumIodidePlusTl->AddMaterial(SODIUM_IODIDE, 99.9*perCent);
	g4SodiumIodidePlusTl->AddMaterial(tl, 0.1*perCent);
	myMaterials.push_back(g4SodiumIodidePlusTl);

	// Deuterated Scintillator from Joey
	//deuterium
	G4Isotope* de = new G4Isotope("de", 1, 2, a = 2.02*g/mole);
	G4Element* d = new G4Element(name="Deuterium",symbol="D",nComponents=1);
	d->AddIsotope(de,fractionMass=100.*perCent);

	//DeuScin
	G4Material* deuScin = new G4Material("BC537", density = 0.954*g/cm3, 3);
	deuScin->AddElement(elH,fractionMass=0.0625*perCent);
	deuScin->AddElement(elC,fractionMass=85.7326*perCent);
	deuScin->AddElement(d,fractionMass=14.2049*perCent);
	myMaterials.push_back(deuScin);

	G4Material* peek = new G4Material("Peek",density = 1.26*g/cm3, nComponents=3);
	peek->AddElement(elC, nAtoms=19);
	peek->AddElement(elH, nAtoms=12);
	peek->AddElement(elO, nAtoms=3);
	myMaterials.push_back(peek);

	G4Material* wTa = new G4Material("WTa",density = 18677*kg/m3, nComponents=2 );
	wTa->AddElement(elW, 80.*perCent);
	wTa->AddElement(elTa, 20.*perCent);
	myMaterials.push_back(wTa);

	G4Material* kapton = new G4Material("Kapton",density = 1.43*g/cm3, nComponents=4);
	kapton->AddElement(elC,nAtoms=22);
	kapton->AddElement(elH,nAtoms=10);
	kapton->AddElement(elN,nAtoms=2);
	kapton->AddElement(elO,nAtoms=5);
	myMaterials.push_back(kapton);

	G4Material* calcium = new G4Material("Calcium", z=20., a=40.078*g/mole, density = 1.55*g/cm3);
	myMaterials.push_back(calcium);

	G4Material* acrylic = new G4Material(name="Acrylic", density=1.18*g/cm3, nElements=3);
	acrylic->AddElement(elC, nAtoms=5);
	acrylic->AddElement(elO, nAtoms=2);
	acrylic->AddElement(elH, nAtoms=8);
	myMaterials.push_back(acrylic);

	G4Material* barium = new G4Material("Barium", z=56., a=137.3*g/mole, density = 3.62*g/cm3);
	myMaterials.push_back(barium);

	G4Material* bismuth = new G4Material("Bismuth", z=83., a=207*g/mole, density = 9.78*g/cm3);
	myMaterials.push_back(bismuth);

	G4Element* elPdO = new G4Element(name="PalladiumOneTen", symbol="Pd", z=46., a=109.91*g/mole);
	myElements.push_back(elPdO);


	G4Material* myPalladium = new G4Material(name="Palladium", density=1.23*g/cm3, nComponents=2);
	myPalladium->AddElement(elPdO, 97.61 *perCent);
	myPalladium->AddElement(elPd, 2.39 *perCent);
	myMaterials.push_back(myPalladium);

	G4Material* Sm = new G4Material(name="Samarium", z=62., a= 152.0*g/mole, density= 7.54*g/cm3);
	myMaterials.push_back(Sm);

    //epoxy-resin , C11H12O3, from https://jlabsvn.jlab.org/svnroot/halla/groups/g2p/HRSMC/src/HRSMaterial.cc
	//density = 0.95*g/cm3;     //volumn ratio 60%:40%
	density = 1.268*g/cm3;  //weight ratio 60%:40%
	G4Material * epoxy = new G4Material(name="Epoxy-Resin", density, nComponents=3);
	epoxy->AddElement(elH, nAtoms=12);
	epoxy->AddElement(elC, nAtoms=11);
	epoxy->AddElement(elO, nAtoms=3);
    myMaterials.push_back(epoxy);

    // make the helium gas material for TI-STAR vacuum chamber
    G4Material* helium = man->FindOrBuildMaterial("G4_He");
    myMaterials.push_back(helium);
    
    // Make the deuterium gas material for TI-STAR target
    G4Material* deuteriumGas = new G4Material("2H", 0.1645*kg/m3, 1, kStateGas, 298*kelvin, 1.*atmosphere);//10.0e-3. * CLHEP::atmosphere, 298 * CLHEP::kelvin
    G4Isotope * iso_H2 = new G4Isotope("H2", 1, 2, 2.014*CLHEP::g/CLHEP::mole);
    G4Element * elD = new G4Element("Deuterium Atom","D",1);
    elD->AddIsotope(iso_H2, 100.*CLHEP::perCent);
    deuteriumGas->AddElement(elD,2);
    myMaterials.push_back(deuteriumGas);
    
    G4Material* xenonGas = new G4Material("Xe", 54, 131.29*g/mole, 5.458*mg/cm3, kStateGas, 293.15*kelvin, 1.*atmosphere);
    myMaterials.push_back(xenonGas);
    
    //G4cout<<*(G4Material::GetMaterialTable())<<G4endl;
}

void DetectorConstruction::DefineTistarTargetMaterials() 
{
	std::vector<G4Material*> myMaterials;   // warnings of unused variables
    
    G4String name = "TistarGasTarget_" + TistarSettings::Get()->GetTargetMaterialName();
    // from TRexMaterials for the TISTAR target
    G4Material * target_material = G4NistManager::Instance()->ConstructNewGasMaterial(name, // material name 
                                                                                      TistarSettings::Get()->GetTargetMaterialName(), // G4Material that this is a gas of - this is why we need to create the 2H G4Material above
                                                                                      TistarSettings::Get()->GetTargetTemperature(),  // temperature
                                                                                      TistarSettings::Get()->GetTargetPressure());    // pressure
    myMaterials.push_back(target_material);
}

void DetectorConstruction::DefineTistarVacuumChamberMaterials()
{
	std::vector<G4Material*> myMaterials;   // warnings of unused variables

    G4String name = "TistarVacuumChamber_" + TistarSettings::Get()->GetVacuumChamberMaterialName();
    G4Material * helium_gas_material = G4NistManager::Instance()->ConstructNewGasMaterial(name, TistarSettings::Get()->GetVacuumChamberMaterialName(), TistarSettings::Get()->GetVacuumChamberGasTemperature(), TistarSettings::Get()->GetVacuumChamberGasPressure());
    myMaterials.push_back(helium_gas_material);
}
