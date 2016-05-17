#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemTestcan.hh"

#include "G4SystemOfUnits.hh"

#include <string>

DetectionSystemTestcan::DetectionSystemTestcan(G4double length, G4double radius) :
    // LogicalVolumes
    testcan_alum_casing_log(0),
    testcan_scintillator_log(0),
    testcan_quartz_window_log(0)
{
    // can properties
    scintillator_length         = length;
    alum_can_thickness          = 1.0*mm;
    scintillator_inner_radius   = 0.0*mm;
    scintillator_outer_radius   = radius;
    quartz_thickness            = 6.35*mm;
    quartz_radius               = radius + alum_can_thickness;
    can_material                = "G4_Al";
    liquid_material             = "Deuterated Scintillator";
    quartz_material             = "G4_SILICON_DIOXIDE";

    start_phi               = 0.0*deg;
    end_phi                 = 360.0*deg;
   
    liquid_colour           = G4Colour(0.0/255.0,255.0/255.0,225.0/255.0);
    grey_colour             = G4Colour(0.5, 0.5, 0.5); 
    quartz_colour           = G4Colour(1.0, 0.0, 1.0);    

}


DetectionSystemTestcan::~DetectionSystemTestcan()
{
    // LogicalVolumes
    delete testcan_alum_casing_log;
    delete testcan_scintillator_log;
    delete testcan_quartz_window_log;
}


G4int DetectionSystemTestcan::Build()
{
    G4AssemblyVolume* myAssemblyTestcan = new G4AssemblyVolume();
    this->assemblyTestcan = myAssemblyTestcan;

    BuildTestcan();

    return 1;
}

G4int DetectionSystemTestcan::PlaceDetector(G4LogicalVolume* exp_hall_log)
{
    G4RotationMatrix * rotate = new G4RotationMatrix;
    G4ThreeVector move = G4ThreeVector(0., 0., 0.);

    assemblyTestcan->MakeImprint(exp_hall_log, move, rotate);

    return 1;
}

G4int DetectionSystemTestcan::BuildTestcan()
{
    G4ThreeVector move, direction;
    G4RotationMatrix* rotate;

    G4Material* can_g4material = G4Material::GetMaterial(this->can_material);
    if( !can_g4material ) {
        G4cout << " ----> Material " << this->can_material << " not found, cannot build! " << G4endl;
        return 0;
    }

    G4Material* liquid_g4material = G4Material::GetMaterial(this->liquid_material);
    if( !liquid_g4material ) {
        G4cout << " ----> Material " << this->liquid_material << " not found, cannot build! " << G4endl;
        return 0;
    }

    G4Material* quartz_g4material = G4Material::GetMaterial(this->quartz_material);
    if(!quartz_g4material ) {
        G4cout << " ----> Material " << this->quartz_material << " not found, cannot build! " << G4endl;
        return 0;
    }


    // SCINTILLATION PROPERTIES
    //-------------------------------------------------------------------------------------------------------------------------
    G4MaterialPropertiesTable * MPT = new G4MaterialPropertiesTable();
    const G4int NUM = 6;
    G4double photon_energies[NUM] = {3.1*CLHEP::eV,2.88*CLHEP::eV,2.82*CLHEP::eV,2.695*CLHEP::eV,2.58*CLHEP::eV,2.48*CLHEP::eV};
    G4double emission_spectra[NUM] = {0.05,1.,0.7,0.37,0.2,0.1};
    MPT->AddConstProperty("RINDEX",1.498);
    MPT->AddConstProperty("SCINTILLATIONYIELD",9200./CLHEP::MeV);

    G4double electron_energy[4] = {1.*CLHEP::keV,1.*CLHEP::MeV,10.*CLHEP::MeV,100*CLHEP::MeV};
    G4double electron[4] = {1.,9200.,92000.,920000.};
     MPT->AddProperty("ELECTRONSCINTILLATIONYIELD",electron_energy,electron,4);
     // MPT->AddProperty("DEUTERONSCINTILLATIONYIELD",1.*CLHEP::MeV,2000./CLHEP::MeV,1);
    /*
    MPT->AddConstProperty("IONSCINTILLATIONYIELD",1000./CLHEP::MeV);
    MPT->AddConstProperty("DEUTERONSCINTILLATIONYIELD",9200./CLHEP::MeV);
    MPT->AddConstProperty("TRITONSCINTILLATIONYIELD",1000./CLHEP::MeV);
    MPT->AddConstProperty("ALPHASCINTILLATIONYIELD",1000./CLHEP::MeV);
    MPT->AddConstProperty("PROTONSCINTILLATIONYIELD",9200./CLHEP::MeV);
    MPT->AddConstProperty("ELECTRONSCINTILLATIONYIELD",1000./CLHEP::MeV);
    */
    MPT->AddConstProperty("RESOLUTIONSCALE",25.);
    MPT->AddConstProperty("ABSLENGTH",3.*CLHEP::m);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // values for fast and slow time component taken from :                                              //
    //     http://research.physics.lsa.umich.edu/twinsol/Publications/Ojaruegathesis.pdf                 //
    //                                                                                                   //
    //     Pulse shape analysis of liquid scintillators for neutron studies,  S. Marrone et al.          //
    MPT->AddConstProperty("FASTTIMECONSTANT",3.5*CLHEP::ns);
    MPT->AddProperty("FASTCOMPONENT",photon_energies,emission_spectra,NUM)->SetSpline(true);
    // MPT->AddConstProperty("SLOWTIMECONSTANT",50.7*CLHEP::ns);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    MPT->AddConstProperty("YIELDRATIO",1.);
    //------------------------------------------------------------------------------------------------------------------------

    liquid_g4material->SetMaterialPropertiesTable(MPT);

    // building the aluminum can
    G4double cut_extra = 10.0*cm;
    G4Tubs* cut = new G4Tubs("cutting volume", this->scintillator_inner_radius, this->scintillator_outer_radius, this->scintillator_length/2.0 + cut_extra, this->start_phi, this->end_phi);
    G4Tubs* can = new G4Tubs("can volume before cut", this->scintillator_inner_radius, this->scintillator_outer_radius + this->alum_can_thickness, this->scintillator_length/2.0 + this->alum_can_thickness/2.0, this->start_phi, this->end_phi);
    
    move = G4ThreeVector(0., 0., -0.5*this->alum_can_thickness - cut_extra);
    rotate = new G4RotationMatrix;
    G4SubtractionSolid* testcan_alum_casing = new G4SubtractionSolid("testcan_alum_can", can, cut, rotate, move);

    // building the scintillator
    G4Tubs* testcan_scintillator = new G4Tubs("testcan_scintillator", this->scintillator_inner_radius, this->scintillator_outer_radius, this->scintillator_length/2.0, this->start_phi, this->end_phi);

    // building the quartz window
    G4Tubs* testcan_quartz_window = new G4Tubs("testcan_quartz_window", this->scintillator_inner_radius, this->quartz_radius, this->quartz_thickness/2.0, this->start_phi, this->end_phi);

    // Set visualization attributes
    G4VisAttributes* can_vis_att = new G4VisAttributes(grey_colour);
    can_vis_att->SetVisibility(true);
    G4VisAttributes* liquid_vis_att = new G4VisAttributes(liquid_colour);
    liquid_vis_att->SetVisibility(true);
    G4VisAttributes* quartz_vis_att = new G4VisAttributes(quartz_colour);
    quartz_vis_att->SetVisibility(true);
    //quartz_vis_att->SetForceWireframe(true);
    
    // Define rotation and movement objects for aluminum can
    direction 	  = G4ThreeVector(0., 0., 1.);
    move          = G4ThreeVector(0., 0., -1.0*this->scintillator_length/2.0 + this->alum_can_thickness/2.0 );
    rotate = new G4RotationMatrix;
    
    //logical volume for aluminum can
    if( testcan_alum_casing_log == NULL )
    {
        testcan_alum_casing_log = new G4LogicalVolume(testcan_alum_casing, can_g4material, "testcan_alum_casing_log", 0, 0, 0);
        testcan_alum_casing_log->SetVisAttributes(can_vis_att);
    }
    this->assemblyTestcan->AddPlacedVolume(testcan_alum_casing_log, move, rotate);
    
    // Define rotation and movement objects for scintillator
    direction 	  = G4ThreeVector(0., 0., 1.);
    move          = G4ThreeVector(0., 0., -1.0*scintillator_length/2.0);
    rotate = new G4RotationMatrix;
    
    // logical volume for scintillator
    if( testcan_scintillator_log == NULL )
    {
        testcan_scintillator_log = new G4LogicalVolume(testcan_scintillator, liquid_g4material, "testcan_scintillator_log", 0, 0, 0);
        testcan_scintillator_log->SetVisAttributes(liquid_vis_att);
    }
    this->assemblyTestcan->AddPlacedVolume(testcan_scintillator_log, move, rotate);

    // Define rotation and movement objects for quartz_window
    direction 	  = G4ThreeVector(0., 0., 1.);
    move          = G4ThreeVector(0., 0., -1.0*quartz_thickness/2.0 - 1.0*scintillator_length);
    rotate = new G4RotationMatrix;
    
    // logical volume for quartz window
    if( testcan_quartz_window_log == NULL )
    {
        testcan_quartz_window_log = new G4LogicalVolume(testcan_quartz_window, quartz_g4material, "testcan_quartz_window_log", 0, 0, 0);
        testcan_quartz_window_log->SetVisAttributes(quartz_vis_att);
    }
    this->assemblyTestcan->AddPlacedVolume(testcan_quartz_window_log, move, rotate);

    
    return 1;
}
