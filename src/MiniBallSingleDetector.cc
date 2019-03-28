/************************************************************************
 * \file MiniBallSingleDetector.hh
 *
 * \class MiniBallSingleDetector   
 * \brief Design of an encapsulated MiniBall HPGe detector
 *
 * \author hans.boie@mpi-hd.mpg.de
 * \author v.bildstein@mpi-hd.mpg.de
 *
 ************************************************************************/

#include "MiniBallSingleDetector.hh"
#include "MiniBallMaterial.hh"

MiniBallSingleDetector::MiniBallSingleDetector(G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(0), phi(0), theta(0),
   spin(0), HPGe_transformation(NULL), capsule(p_capsule),
   histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(G4String p_name,
  G4VPhysicalVolume* p_mother, G4bool p_capsule) 
  :name(p_name), mother(p_mother), rho(0), phi(0), theta(0),
   spin(0), HPGe_transformation(NULL), capsule(p_capsule),
   histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(G4String p_name,
  G4VPhysicalVolume* p_mother, G4double p_rho, G4double p_theta,
  G4double p_phi, G4bool p_capsule)
  :name(p_name), mother(p_mother), rho(p_rho), phi(p_phi),
   theta(p_theta), spin(0), HPGe_transformation(NULL), capsule(p_capsule),
   histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(G4String p_name,
  G4VPhysicalVolume* p_mother, G4double p_rho, G4double p_theta,
  G4double p_phi, G4double p_spin, G4bool p_capsule)
  :name(p_name), mother(p_mother), rho(p_rho), phi(p_phi),
   theta(p_theta), spin(p_spin), HPGe_transformation(NULL), 
   capsule(p_capsule), histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(G4String p_name,
  G4VPhysicalVolume* p_mother, G4Transform3D* p_transform, G4bool p_capsule) 
  :name(p_name), mother(p_mother), rho(0), phi(0), theta(0),
   spin(0), HPGe_transformation(p_transform), capsule(p_capsule),
   histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(G4String p_name,
  G4VPhysicalVolume* p_mother, G4Transform3D* p_transform,
  G4double p_spin, G4bool p_capsule) 
  :name(p_name), mother(p_mother), rho(0), phi(0), theta(0), 
   spin(p_spin), HPGe_transformation(p_transform), capsule(p_capsule),
   histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(G4double p_rho,
  G4double p_theta, G4double p_phi, G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(p_rho), phi(p_phi),
   theta(p_theta), spin(0), HPGe_transformation(NULL), capsule(p_capsule),
   histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(G4double p_rho, 
  G4double p_theta, G4double p_phi, G4double p_spin, G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(p_rho), phi(p_phi),
   theta(p_theta), spin(p_spin), HPGe_transformation(NULL), 
   capsule(p_capsule), histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(G4Transform3D* p_transform, 
    G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(0), phi(0), theta(0),
   spin(0), HPGe_transformation(p_transform), capsule(p_capsule),
   histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(G4Transform3D* p_transform,
  G4double p_spin, G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(0), phi(0), theta(0),
   spin(p_spin), HPGe_transformation(p_transform), capsule(p_capsule),
   histoManager(0),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(0), phi(0), theta(0),
   spin(0), HPGe_transformation(NULL), capsule(p_capsule),
   histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4String p_name, G4VPhysicalVolume* p_mother, G4bool p_capsule) 
  :name(p_name), mother(p_mother), rho(0), phi(0),
   theta(0), spin(0), HPGe_transformation(NULL), capsule(p_capsule),
   histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4String p_name,G4VPhysicalVolume* p_mother, G4double p_rho, 
  G4double p_theta, G4double p_phi, G4bool p_capsule)
  :name(p_name), mother(p_mother), rho(p_rho), phi(p_phi),
   theta(p_theta), spin(0), HPGe_transformation(NULL), capsule(p_capsule),
   histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4String p_name, G4VPhysicalVolume* p_mother, G4double p_rho, 
  G4double p_theta, G4double p_phi, G4double p_spin, G4bool p_capsule)
  :name(p_name), mother(p_mother), rho(p_rho), phi(p_phi),
   theta(p_theta), spin(p_spin), HPGe_transformation(NULL), 
   capsule(p_capsule), histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4String p_name, G4VPhysicalVolume* p_mother, G4Transform3D* p_transform,
  G4bool p_capsule) 
  :name(p_name), mother(p_mother), rho(0), phi(0),
   theta(0), spin(0), HPGe_transformation(p_transform), capsule(p_capsule),
   histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4String p_name, G4VPhysicalVolume* p_mother, G4Transform3D* p_transform,
  G4double p_spin, G4bool p_capsule) 
  :name(p_name), mother(p_mother), rho(0), phi(0),
   theta(0), spin(p_spin), HPGe_transformation(p_transform), 
   capsule(p_capsule), histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE),
   det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4double p_rho, G4double p_theta, G4double p_phi, G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(p_rho), phi(p_phi),
   theta(p_theta), spin(0), HPGe_transformation(NULL), capsule(p_capsule),
   histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE), det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4double p_rho, G4double p_theta,G4double p_phi, G4double p_spin, 
  G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(p_rho), phi(p_phi),
   theta(p_theta), spin(p_spin), HPGe_transformation(NULL), 
   capsule(p_capsule), histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE), det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4Transform3D* p_transform, G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(0), phi(0), theta(0),
   spin(0), HPGe_transformation(p_transform), capsule(p_capsule),
   histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE), det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::MiniBallSingleDetector(MiniBallHistoManager* HM,
  G4Transform3D* p_transform, G4double p_spin, G4bool p_capsule)
  :name("/HPGe"), mother(NULL), rho(0), phi(0), theta(0),
   spin(p_spin), HPGe_transformation(p_transform), capsule(p_capsule),
   histoManager(HM),  eff_fact(1.),
   deadlayer_thickness(DETECTOR_DEFAULT_DEADLAYER_THICKNESS),
   deadlayer_type(DEFAULT_DEADLAYER_TYPE), det_vis_att(NULL), capsule_vis_att(NULL), 
   capsule_cut(false), det_cut(false)
{}

MiniBallSingleDetector::~MiniBallSingleDetector() 
{; }

// construction of the detector, if no mother volume is given an world
// volume is created, gives back either the mothervolume with the
// detector is placed in or the new created world volume

G4VPhysicalVolume* MiniBallSingleDetector::Construct()
{
  if(histoManager)
    {
      G4int det_nb = histoManager->GetNbOfDetectors();

		std::ostringstream tmp_strm;
      tmp_strm.clear();
      tmp_strm << det_nb;
      name += "_" + tmp_strm.str();
    }

  if(name.at(0) != '/')
    name = "/"+name;

  // germanium
  G4Material* Ge = MiniBallMaterial::Get()->GetMaterialGe();
  
  // aluminum
  G4Material* Al = MiniBallMaterial::Get()->GetMaterialAl();

  G4Material* vacuum = MiniBallMaterial::Get()->GetMaterialVacuum();

#ifdef G4VIS_USE  
  // visualisation attributes
  if(!det_vis_att)
    {
      G4VisAttributes* HPGeVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
      HPGeVisAtt->SetVisibility(true);
      det_vis_att = HPGeVisAtt;
    }
  if(!capsule_vis_att)
    {
      G4VisAttributes* AlVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
      AlVisAtt->SetVisibility(true);
      capsule_vis_att = AlVisAtt;
    }
#endif

  // construct germanium crystal
  G4double det_visangle=360.*deg;
  G4double det_startangle=0.*deg;
  if(det_cut)
    {
      det_visangle=180.*deg;
      if(det_cut_type==0)
	det_startangle=180.*deg;
    }

  G4double crystal_height = DETECTOR_HEIGHT*eff_fact;
  G4double crystal_min_radius = DETECTOR_MIN_RADIUS*eff_fact;
  G4double crystal_max_radius = DETECTOR_MAX_RADIUS*eff_fact;

  G4LogicalVolume* HPGe_log;
  G4LogicalVolume* deadlayer_log;

  if(deadlayer_type == DEADLAYER_FRONT)
    {
      G4double z[9] =    {crystal_height/2-deadlayer_thickness,
			  24.*mm, 23.9*mm, 23.5*mm, 23.*mm,
			  22.*mm, 21.*mm, 20.*mm, -crystal_height/2};
      G4double rmin[9] = { 0.*mm,  0.*mm,   1.*mm,   2.*mm,  3.*mm,   
			   4.5*mm, 4.9*mm,  5.*mm,   5.*mm};
      G4double rmax[9] = {crystal_max_radius, crystal_max_radius,
			  crystal_max_radius, crystal_max_radius,
			  crystal_max_radius, crystal_max_radius,
			  crystal_max_radius, crystal_max_radius,
			  crystal_max_radius};

      G4Polycone* detector_tub;
      if(!det_cut)
	detector_tub
	  = new G4Polycone(name+"/detector_tub", 0.*deg, 360.*deg, 9,
			   z, rmin, rmax);
      else
	detector_tub
	  = new G4Polycone(name+"/detector_tub", det_startangle, det_visangle,
			   9, z, rmin, rmax);

      // cutting the edge
      G4double c_z[] =    {crystal_height/2, -crystal_height/2};
      
      G4double c_rmin[] = {0.*mm, 0.*mm};
  
      G4double c_rmax[] = {crystal_min_radius, crystal_max_radius};
  
      G4Polyhedra* detector_cut
	= new G4Polyhedra(name+"/detector_edge", 0.*deg, 360.*deg, 6, 2,
			  c_z, c_rmin, c_rmax);
      
      // tub with edges cutted off
      G4IntersectionSolid* HPGe_sol
	= new G4IntersectionSolid(name+"/HPGe_sol", detector_tub, detector_cut,
				  0, G4ThreeVector());
  
      HPGe_log = new G4LogicalVolume(HPGe_sol, Ge, name+"/HPGe_log");
      
#ifdef G4VIS_USE  
      HPGe_log->SetVisAttributes(det_vis_att);
#endif
    
      //dead layer for the detector
      G4double dead_z[] = {crystal_height/2,
			   crystal_height/2-deadlayer_thickness};
      G4double dead_rmin[] = {0, 0};
      G4double dead_rmax[] = {crystal_min_radius, crystal_min_radius};
      
      G4Polyhedra* deadlayer_sol;
      
      if(!det_cut)
	deadlayer_sol
	  = new G4Polyhedra(name+"/deadlayer_sol", 0*deg, 360*deg, 6, 2,
			    dead_z, dead_rmin, dead_rmax);
      else
	deadlayer_sol
	  = new G4Polyhedra(name+"/deadlayer_sol", det_startangle,
			    det_visangle, 3, 2,
			    dead_z, dead_rmin, dead_rmax);
      
      deadlayer_log
	= new G4LogicalVolume(deadlayer_sol, Ge, name+"/deadlayer_log");
      
#ifdef G4VIS_USE
      deadlayer_log->SetVisAttributes(det_vis_att);
#endif
    }
  else if(deadlayer_type == DEADLAYER_OUTER_CONTACT)
    {
      G4double z[9] =    {crystal_height/2-deadlayer_thickness,
			  24.*mm, 23.9*mm, 23.5*mm, 23.*mm,
			  22.*mm, 21.*mm, 20.*mm, -crystal_height/2};
      G4double rmin[9] = { 0.*mm,  0.*mm,   1.*mm,   2.*mm,  3.*mm,   
			   4.5*mm, 4.9*mm,  5.*mm,   5.*mm};
      G4double rmax[9] = {crystal_max_radius-deadlayer_thickness,
			  crystal_max_radius-deadlayer_thickness,
			  crystal_max_radius-deadlayer_thickness,
			  crystal_max_radius-deadlayer_thickness,
			  crystal_max_radius-deadlayer_thickness,
			  crystal_max_radius-deadlayer_thickness,
			  crystal_max_radius-deadlayer_thickness,
			  crystal_max_radius-deadlayer_thickness,
			  crystal_max_radius-deadlayer_thickness};

      G4Polycone* detector_tub;
      if(!det_cut)
	detector_tub
	  = new G4Polycone(name+"/detector_tub", 0.*deg, 360.*deg, 9,
			   z, rmin, rmax);
      else
	detector_tub
	  = new G4Polycone(name+"/detector_tub", det_startangle, det_visangle,
			   9, z, rmin, rmax);

      // cutting the edge
      G4double c_z[] =    {crystal_height/2, -crystal_height/2};
      
      G4double c_rmin[] = {0.*mm, 0.*mm};
  
      G4double c_rmax[] = {crystal_min_radius-deadlayer_thickness,
			   crystal_max_radius-deadlayer_thickness};
  
      G4Polyhedra* detector_cut
	= new G4Polyhedra(name+"/detector_edge", 0.*deg, 360.*deg, 6, 2,
			  c_z, c_rmin, c_rmax);
      
      // tub with edges cutted off
      G4IntersectionSolid* HPGe_sol
	= new G4IntersectionSolid(name+"/HPGe_sol", detector_tub, detector_cut,
				  0, G4ThreeVector());
  
      HPGe_log = new G4LogicalVolume(HPGe_sol, Ge, name+"/HPGe_log");
      
#ifdef G4VIS_USE  
      HPGe_log->SetVisAttributes(det_vis_att);
#endif

      //dead layer for the detector
      G4double dead_z[2] =    {crystal_height/2,
			       -crystal_height/2};
      G4double dead_rmin[2] = {0.*mm, 0.*mm};
      G4double dead_rmax[2] = {crystal_max_radius, crystal_max_radius};

      G4Polycone* deadlayer_tub;
      if(!det_cut)
	deadlayer_tub
	  = new G4Polycone(name+"/deadlayer_tub", 0.*deg, 360.*deg, 2,
			   dead_z, dead_rmin, dead_rmax);
      else
	deadlayer_tub
	  = new G4Polycone(name+"/deadlayer_tub", det_startangle, det_visangle,
			   2, dead_z, dead_rmin, dead_rmax);

      // cutting the edge
      c_rmax[0] = crystal_min_radius;
      c_rmax[1] = crystal_max_radius;
  
      G4Polyhedra* deadlayer_cut
	= new G4Polyhedra(name+"/deadlayer_edge", 0.*deg, 360.*deg, 6, 2,
			  c_z, c_rmin, c_rmax);
      
      // tub with edges cutted off
      G4IntersectionSolid* deadlayer_cut_sol
	= new G4IntersectionSolid(name+"/deadlayer_cut_sol", deadlayer_tub,
				  deadlayer_cut,
				  0, G4ThreeVector());

      // cut out detector volume
      dead_z[0] = crystal_height/2-deadlayer_thickness;
      dead_z[1] = -crystal_height/2-SMALL_DIST;
      dead_rmax[0] = crystal_max_radius-deadlayer_thickness;
      dead_rmax[1] = crystal_max_radius-deadlayer_thickness;

      G4Polycone* deadlayer_cut_tub
	= new G4Polycone(name+"/deadlayer_cut_tub", 0.*deg, 360.*deg, 2,
			 dead_z, dead_rmin, dead_rmax);

      // cutting the edge
      c_rmax[0] = crystal_min_radius-deadlayer_thickness;
      c_rmax[1] = crystal_max_radius-deadlayer_thickness;
  
      G4Polyhedra* deadlayer_cut_edge
	= new G4Polyhedra(name+"/deadlayer_cut_dege", 0.*deg, 360.*deg, 6, 2,
			  c_z, c_rmin, c_rmax);

      // placeholder for detector volume
      G4IntersectionSolid* deadlayer_cut_sol2
	= new G4IntersectionSolid(name+"/deadlayer_cut_sol2",
				  deadlayer_cut_tub, deadlayer_cut_edge);
      
      G4SubtractionSolid* deadlayer_sol
	= new G4SubtractionSolid(name+"/deadlayer_sol", deadlayer_cut_sol,
	                         deadlayer_cut_sol2,
				 0, G4ThreeVector());
      deadlayer_log
	= new G4LogicalVolume(deadlayer_sol, Ge, name+"/deadlayer_log");
      
#ifdef G4VIS_USE
      deadlayer_log->SetVisAttributes(G4VisAttributes::Invisible);
#endif
    }
  else
    {
      G4cerr << "unknown deadlayer_type: " << deadlayer_type
	     << ", exiting..." << G4endl;
      exit(-1);
    }

  // construct cap of capsule
  std::vector<G4LogicalVolume*> slice_frontA_log;
  std::vector<G4Transform3D> slice_frontA_transform;
  std::vector<G4LogicalVolume*> slice_frontB_log;
  std::vector<G4Transform3D> slice_frontB_transform;
  std::vector<G4LogicalVolume*> slice_back_log;
  std::vector<G4Transform3D> slice_back_transform;

  G4LogicalVolume* capsule_front_log = NULL;
  G4LogicalVolume* capsule_back_log = NULL;
   
  if(capsule)
    {
      G4double cap_visangle=360.*deg;
      G4double cap_startangle=0.*deg;
      if(capsule_cut)
	{
	  cap_visangle=180.*deg;
	  if(capsule_cut_type==0)
	    cap_startangle=180.*deg;
	}
      
      // construct the front cap
      G4double front_z[] = {CAPSULE_FRONT_THICKNESS/2,
			    -CAPSULE_FRONT_THICKNESS/2};
      
      G4double front_rmin[] = {0.*mm, 0.*mm};
      
      G4double front_rmax[] = {SW_DETECTOR_FRONT/2
			       +(GAP_TO_CAPSULE+CAPSULE_THICKNESS)
			       *(1/cos(DETECTOR_ANGLE)-tan(DETECTOR_ANGLE)),
			       SW_DETECTOR_FRONT/2
			       +(GAP_TO_CAPSULE+CAPSULE_THICKNESS)
			       /cos(DETECTOR_ANGLE)
			       -GAP_TO_CAPSULE*tan(DETECTOR_ANGLE)};
      
      G4Polyhedra* capsule_front;
      if(!capsule_cut)
	capsule_front
	  = new G4Polyhedra(name+"/capsule_front", 0.*deg, 360.*deg, 6, 2,
			    front_z, front_rmin, front_rmax);
      else
	capsule_front
	  = new G4Polyhedra(name+"/capsule_front",
			    cap_startangle, cap_visangle,
			    3, 2, front_z, front_rmin, front_rmax);
      
      capsule_front_log
	= new G4LogicalVolume(capsule_front, Al, name+"/capsule_front_log");
      
#ifdef G4VIS_USE
      capsule_front_log->SetVisAttributes(capsule_vis_att);
#endif
      
      // construct back of cap
      
      G4double back_z[] = {3*CAPSULE_THICKNESS/2+CAPSULE_CONTACT_HOLE_DEPTH,
			   CAPSULE_THICKNESS/2+CAPSULE_CONTACT_HOLE_DEPTH,
			   CAPSULE_THICKNESS/2+CAPSULE_CONTACT_HOLE_DEPTH,
			   CAPSULE_THICKNESS/2,
			   CAPSULE_THICKNESS/2,
			   -CAPSULE_THICKNESS/2};
      
      G4double back_rmin[] = {0.*mm,
			      0.*mm,
			      CAPSULE_CONTACT_HOLE_RAD,
			      CAPSULE_CONTACT_HOLE_RAD,
			      CAPSULE_CONTACT_HOLE_RAD,
			      CAPSULE_CONTACT_HOLE_RAD};
      
      G4double back_rmax[] = {CAPSULE_CONTACT_HOLE_RAD+CAPSULE_THICKNESS,
			      CAPSULE_CONTACT_HOLE_RAD+CAPSULE_THICKNESS,
			      CAPSULE_CONTACT_HOLE_RAD+CAPSULE_THICKNESS,
			      CAPSULE_CONTACT_HOLE_RAD+CAPSULE_THICKNESS,
			      DETECTOR_MAX_RADIUS+GAP_TO_CAPSULE,
			      DETECTOR_MAX_RADIUS+GAP_TO_CAPSULE};
      
      G4Polyhedra* capsule_back;
      if(!capsule_cut)
	capsule_back
	  = new G4Polyhedra(name+"/capsule_back", 0.*deg, 360.*deg,
			    NB_OF_SLICES*6, 6,
			    back_z, back_rmin, back_rmax);
      else
	capsule_back
	  = new G4Polyhedra(name+"/capsule_back", cap_startangle, cap_visangle,
			    NB_OF_SLICES*3, 6,
			    back_z, back_rmin, back_rmax);
      
      capsule_back_log
	= new G4LogicalVolume(capsule_back, Al, name+"/capsule_back_log");
      
#ifdef G4VIS_USE
      capsule_back_log->SetVisAttributes(capsule_vis_att);
#endif
      // build up the rest of the cap out of slices
      G4double slice_begin;
      G4double slice_end = -30.*deg;
      
      G4RotationMatrix rot_slice;
      rot_slice.rotateY(M_PI);
      rot_slice.rotateZ(30.*deg);
      G4Transform3D spin_slice(rot_slice, G4ThreeVector(0., 0., 0.));
      
      for(int slice=0; slice < NB_OF_SLICES; slice++)
	{
	  int side;
	  
	  slice_begin = slice_end;
	  slice_end += 60.*deg/NB_OF_SLICES;
	  
	  // definition of the intersecting planes
	  
	  MiniBallPlane begin_plane(G4ThreeVector(0, 0, 0),
				    G4ThreeVector(0, 0, 1),
				    G4ThreeVector(1, tan(slice_begin), 0));
	  
	  MiniBallPlane end_plane(G4ThreeVector(0, 0, 0),
				  G4ThreeVector(0, 0, 1),
				  G4ThreeVector(1, tan(slice_end), 0));
	  
	  MiniBallPlane front_cut_plane(G4ThreeVector(
		    DETECTOR_MIN_RADIUS+GAP_TO_CAPSULE*cos(DETECTOR_ANGLE),
                    0,
		   -CAPSULE_HEIGHT/2-GAP_TO_CAPSULE*sin(DETECTOR_ANGLE)),
		    G4ThreeVector(0, 1, 0),
		    G4ThreeVector(tan(DETECTOR_ANGLE), 0, 1));

	  MiniBallPlane front_cut_plane_top(G4ThreeVector(
		    DETECTOR_MIN_RADIUS+
		    (GAP_TO_CAPSULE+CAPSULE_THICKNESS)*cos(DETECTOR_ANGLE),
                    0,
		   -CAPSULE_HEIGHT/2-
		    (GAP_TO_CAPSULE+CAPSULE_THICKNESS)*sin(DETECTOR_ANGLE)),
			        G4ThreeVector(0, 1, 0),
				G4ThreeVector(tan(DETECTOR_ANGLE), 0, 1));

	  G4double slice_mid = (slice_begin+slice_end)/2;
	  G4double back_cut_plane_position_length =
	    DETECTOR_MAX_RADIUS+GAP_TO_CAPSULE;
	  MiniBallPlane back_cut_plane(G4ThreeVector(
		    back_cut_plane_position_length*cos(slice_mid),	
		    back_cut_plane_position_length*sin(slice_mid),	
		    0),
				 G4ThreeVector(0, 0, 1),
				 G4ThreeVector(-sin(slice_mid),	
		                               cos(slice_mid), 0));

	  G4double back_cut_plane_position_length_top =
	    DETECTOR_MAX_RADIUS+GAP_TO_CAPSULE+CAPSULE_THICKNESS;
	  MiniBallPlane back_cut_plane_top(G4ThreeVector(
		    back_cut_plane_position_length_top*cos(slice_mid),	
		    back_cut_plane_position_length_top*sin(slice_mid),	
		    0),
				 G4ThreeVector(0, 0, 1),
				 G4ThreeVector(-sin(slice_mid),	
		                               cos(slice_mid), 0));

	  MiniBallPlane front_plane(G4ThreeVector(0, 0, 
				    -CAPSULE_HEIGHT/2-GAP_TO_CAPSULE),
				    G4ThreeVector(1, 0, 0),
				    G4ThreeVector(0, 1, 0));

	  MiniBallPlane back_plane(G4ThreeVector(0, 0, CAPSULE_HEIGHT/2),
				   G4ThreeVector(1, 0, 0),
				   G4ThreeVector(0, 1, 0));

	  // cutting lines (bottom)
	  MiniBallLine cut_line = front_cut_plane.Intersect(back_cut_plane);
	  MiniBallLine front_line = front_plane.Intersect(front_cut_plane);
	  MiniBallLine back_line = back_plane.Intersect(back_cut_plane);
	   
	  G4ThreeVector x[8];
	  G4ThreeVector y[8];
	  
	  // slices coordinates (bottom)
	  x[0] = begin_plane.Intersect(front_line);
	  x[1] = end_plane.Intersect(front_line);
	  x[2] = end_plane.Intersect(cut_line);
	  x[3] = begin_plane.Intersect(cut_line);

	  y[0] = x[2];
	  y[1] = end_plane.Intersect(back_line);
	  y[2] = x[3];
	  y[3] = begin_plane.Intersect(back_line);

	  if(y[0].z() > y[1].z() || y[2].z() > y[3].z())
	    std::cout << "coordinate error increase NB_OF_SLICES" << std::endl;
	  
	  // cutting lines (top)
	  MiniBallLine cut_line_top =
	    front_cut_plane_top.Intersect(back_cut_plane_top);
	  MiniBallLine front_line_top =
	    front_plane.Intersect(front_cut_plane_top);
	  MiniBallLine back_line_top =
	    back_plane.Intersect(back_cut_plane_top);
	  
	  // slices coordinates (top)
	  x[4] = begin_plane.Intersect(front_line_top);
	  x[5] = end_plane.Intersect(front_line_top);
	  x[6] = end_plane.Intersect(cut_line_top);
	  x[7] = begin_plane.Intersect(cut_line_top);

	  y[4] = x[6];
	  y[5] = end_plane.Intersect(back_line_top);
	  y[6] = x[7];
	  y[7] = begin_plane.Intersect(back_line_top);

	  char slice_name[30];
	  sprintf(slice_name, "/slice_back%d", slice);
	  
	  G4double Dz, Theta, Phi;
	  G4double Dy1, Dx1, Dx2, Alp1;
	  G4double Dy2, Dx3, Dx4, Alp2;

	  G4Transform3D sf_trans
	    = MiniBallGeom::GetTrap(y, &Dz, &Theta, &Phi,
				    &Dy1, &Dx1, &Dx2, &Alp1,
				    &Dy2, &Dx3, &Dx4, &Alp2);

	  slice_back_transform.push_back(spin_slice * sf_trans);

	  G4Trap* new_slice
	    = new G4Trap(name+slice_name, Dz, Theta, Phi,
			 Dy1, Dx1, Dx2, Alp1, Dy2, Dx3, Dx4, Alp2);

	  for(side=1; side<7; side++)
	    {
	      sprintf(slice_name, "/slice_back_log%d_side%d", slice, side);
	      
	      G4LogicalVolume* new_slice_log
		= new G4LogicalVolume(new_slice, Al, name+slice_name);
	      
#ifdef G4VIS_USE
	      new_slice_log->SetVisAttributes(capsule_vis_att);
#endif
	      
	      slice_back_log.push_back(new_slice_log);
	    }
	  
	  G4ThreeVector x1[8];
	  G4ThreeVector x2[8];
	  
	  G4bool trap_two = true;
	  
	  G4ThreeVector tr2_cut_pos;
	  G4ThreeVector tr2_cut_dir;
	  G4ThreeVector tr2_end_pos;
	  G4ThreeVector tr2_end_dir;
	  G4bool tr2_invert = false;
	  
	  if( abs(x[2].z() - x[3].z()) < TRAP_THRESH )
	    {
	      trap_two = false;
	      x1[0] = x[0];
	      x1[1] = x[1];
	      x1[2] = x[3];
	      x1[3] = x[2];
	      x1[4] = x[4];
	      x1[5] = x[5];
	      x1[6] = x[7];
	      x1[7] = x[6];
	    }
	  else if( x[2].z() > x[3].z() )
	    {
	      tr2_cut_pos = x[3];
	      tr2_cut_dir = x[7]-x[3];
	      tr2_end_pos = x[2];
	      tr2_end_dir = x[6]-x[2];
	      tr2_invert = true;
	    }
	  else
	    {
	      tr2_cut_pos = x[2];
	      tr2_cut_dir = x[6]-x[2];
	      tr2_end_pos = x[3];
	      tr2_end_dir = x[7]-x[3];
	      tr2_invert = false;
	    }
	  
	  if(trap_two)
	    {
	      MiniBallPlane tr2_cut_plane(tr2_cut_pos,
					  G4ThreeVector(0, 1, 0),
					  tr2_cut_dir);
	      MiniBallLine tr2_cut_line
		= tr2_cut_plane.Intersect(front_cut_plane);
	      MiniBallLine tr2_cut_line_top
		= tr2_cut_plane.Intersect(front_cut_plane_top);
	      
	      G4ThreeVector c1 = begin_plane.Intersect(tr2_cut_line);
	      G4ThreeVector c2 = end_plane.Intersect(tr2_cut_line);
	      G4ThreeVector c3 = begin_plane.Intersect(tr2_cut_line_top);
	      G4ThreeVector c4 = end_plane.Intersect(tr2_cut_line_top);
	      
	      x1[0] = x[0];
	      x1[1] = x[1];
	      x1[2] = c1;
	      x1[3] = c2;
	      x1[4] = x[4];
	      x1[5] = x[5];
	      x1[6] = c3;
	      x1[7] = c4;
	      
	      tr2_end_pos.setZ(tr2_end_pos.z()-0.001*mm);
	      
	      MiniBallPlane tr2_end_plane(tr2_end_pos,
					  G4ThreeVector(0, 1, 0),
					  tr2_end_dir);
	      MiniBallLine tr2_end_line
		= tr2_end_plane.Intersect(front_cut_plane);
	      MiniBallLine tr2_end_line_top
		= tr2_end_plane.Intersect(front_cut_plane_top);
	      
	      G4ThreeVector e1 = begin_plane.Intersect(tr2_end_line);
	      G4ThreeVector e2 = end_plane.Intersect(tr2_end_line);
	      G4ThreeVector e3 = begin_plane.Intersect(tr2_end_line_top);
	      G4ThreeVector e4 = end_plane.Intersect(tr2_end_line_top);
	      
	      G4ThreeVector e5 = tr2_end_plane.Intersect(cut_line);
	      G4ThreeVector e6 = tr2_end_plane.Intersect(cut_line_top);
	      
	      if(tr2_invert)
		{
		  x2[0] = c1;
		  x2[1] = c2;
		  x2[2] = e5;
		  x2[3] = e2;
		  x2[4] = c3;
		  x2[5] = c4;
		  x2[6] = e6;
		  x2[7] = e4;
		}
	      else
		{
		  x2[0] = c1;
		  x2[1] = c2;
		  x2[2] = e1;
		  x2[3] = e5;
		  x2[4] = c3;
		  x2[5] = c4;
		  x2[6] = e3;
		  x2[7] = e6;
		}		 
	    }
	  
	  sprintf(slice_name, "/slice_frontA_%d", slice);

	  sf_trans
	    = MiniBallGeom::GetTrap(x1, &Dz, &Theta, &Phi,
				    &Dy1, &Dx1, &Dx2, &Alp1,
				    &Dy2, &Dx3, &Dx4, &Alp2);
	  slice_frontA_transform.push_back(spin_slice * sf_trans);

	  new_slice
	    = new G4Trap(name+slice_name, Dz, Theta, Phi,
			 Dy1, Dx1, Dx2, Alp1, Dy2, Dx3, Dx4, Alp2);

	  for(side=1; side<7; side++)
	    {
	      sprintf(slice_name, "/slice_frontA_log%d_side%d", slice, side);
	      
	      G4LogicalVolume* new_slice_log
		= new G4LogicalVolume(new_slice, Al, name+slice_name);
	      
#ifdef G4VIS_USE
	      new_slice_log->SetVisAttributes(capsule_vis_att);
#endif

	      slice_frontA_log.push_back(new_slice_log);
	    }
	  
	  if(trap_two)
	    {
	      sprintf(slice_name, "/slice_frontB_%d", slice);
	      
	      sf_trans
		= MiniBallGeom::GetTrap(x2, &Dz, &Theta, &Phi,
					&Dy1, &Dx1, &Dx2, &Alp1,
					&Dy2, &Dx3, &Dx4, &Alp2);
	      slice_frontB_transform.push_back(spin_slice * sf_trans);
	      
	      new_slice
		= new G4Trap(name+slice_name, Dz, Theta, Phi,
			     Dy1, Dx1, Dx2, Alp1, Dy2, Dx3, Dx4, Alp2);
	      
	      for(side=1; side<7; side++)
		{
		  sprintf(slice_name, "/slice_frontB_log%d_side%d",
			  slice, side);
		  
		  G4LogicalVolume* new_slice_log
		    = new G4LogicalVolume(new_slice, Al, name+slice_name);
		  
#ifdef G4VIS_USE
		  new_slice_log->SetVisAttributes(capsule_vis_att);
#endif
		  
		  slice_frontB_log.push_back(new_slice_log);
		}
	    }
	  else
	    {
	      G4RotationMatrix null_rot;
	      G4Transform3D null_trans(null_rot, G4ThreeVector(0.,0.,0.));
	      slice_frontB_transform.push_back(null_trans);
	      for(side=1; side<7; side++)
		slice_frontB_log.push_back(NULL);
	    }	       
     	}
    }
  
  // contruct mother volume for detector
  G4double env_min_rad = SW_DETECTOR_FRONT/2+
    +(SMALL_DIST+GAP_TO_CAPSULE+CAPSULE_FRONT_THICKNESS)
    *(1/cos(DETECTOR_ANGLE)-tan(DETECTOR_ANGLE));
  
  G4double env_length = GetDetectorEnvelopeLength();
  
  G4double env_max_rad = env_min_rad + env_length*tan(DETECTOR_ANGLE);
  
  G4double env_z[] = {GetDetectorEnvelopeMidPosition(),
		      GetDetectorEnvelopeMidPosition()
		      -GetDetectorEnvelopeLength()};
  
  G4double env_rmin[] = {0.*mm, 0.*mm};
  
  G4double env_rmax[] = {env_min_rad, env_max_rad};
  
  G4Polyhedra* det_env_sol
    = new G4Polyhedra(name+"/det_env", 0.*deg, 360.*deg, 6, 2,
		      env_z, env_rmin, env_rmax);
  
  G4LogicalVolume* det_env_log
    = new G4LogicalVolume(det_env_sol, vacuum, name+"/det_env_log");
  
#ifdef G4VIS_USE
  // mother volume is invisible
  det_env_log->SetVisAttributes(G4VisAttributes::Invisible);
#endif
  
  // if no mothervolume is specified in constructor, create world volume
  if(mother == NULL)
    {
      MiniBallWorld* world = new MiniBallWorld();
      mother = world->Construct();
    }
    
  // set position of detector in mother volume
  if(HPGe_transformation==NULL)
    {
      G4RotationMatrix HPGe_rotation;
      HPGe_rotation.rotateZ(spin);
      HPGe_rotation.rotateX(M_PI);
      HPGe_rotation.rotateY(theta);
      HPGe_rotation.rotateZ(phi);
      
      //distance from the origin to the center of the detector
      G4double trans_mag = rho + MoveToDetectorFront();
      
      G4ThreeVector HPGe_translation(0., 0., trans_mag);
      HPGe_translation.rotateY(theta);
      HPGe_translation.rotateZ(phi);
      
      HPGe_transformation
	= new G4Transform3D(HPGe_rotation, HPGe_translation);
    }
  else if(spin != 0)
    {
      G4RotationMatrix HPGe_rotation = HPGe_transformation->getRotation();
      G4ThreeVector HPGe_translation = HPGe_transformation->getTranslation();
      G4RotationMatrix spin_rotation;
      spin_rotation.rotateZ(spin);
      
      HPGe_rotation *= spin_rotation;
      
      HPGe_transformation 
	= new G4Transform3D(HPGe_rotation, HPGe_translation);
    }
  
  det_env = new G4PVPlacement(*HPGe_transformation, name+"/env",
			      det_env_log, mother, false, 0);
  
  HPGe = new G4PVPlacement(0, G4ThreeVector(0, 0,
				        (CAPSULE_HEIGHT-DETECTOR_HEIGHT)/2),
			   name+"/crystal", HPGe_log, det_env, false, 0);

  deadlayer = new G4PVPlacement(0, G4ThreeVector(0, 0,
				       (CAPSULE_HEIGHT-DETECTOR_HEIGHT)/2),
				name+"/deadlayer",
				deadlayer_log, det_env, false, 0);
 
  if(capsule)
    {
      //set position of capsule front in mother volumen
      G4double trans_mag
	= -CAPSULE_HEIGHT/2-CAPSULE_FRONT_THICKNESS/2-GAP_TO_CAPSULE;
      G4ThreeVector front_translation(0.,0.,-trans_mag);
      
      new G4PVPlacement(0, front_translation, name+"/capsule_front",
			capsule_front_log, det_env, false, 0);
      
      trans_mag
	= -CAPSULE_HEIGHT/2+CAPSULE_BACK;
      G4ThreeVector back_translation(0.,0.,-trans_mag);
      
      new G4PVPlacement(0, back_translation, name+"/capsule_back",
			capsule_back_log, det_env, false, 0);
      
      for(int slice=0; slice < NB_OF_SLICES; slice++)
	{
	  G4int nb_of_sides=6;
	  G4int start_side=0;
	  if(capsule_cut)
	    {
	      if(capsule_cut_type==1)
		start_side=3;
	      else
		nb_of_sides=3;
	    }		  
	  
	  for(int side=start_side; side<nb_of_sides; side++)
	    {
	      char sl_name[255];
	      sprintf(sl_name, "/%s_cap_slice_back_%d_%d",
		      name.c_str(), slice, side);
	      
	      G4RotationMatrix rot_sl;
	      rot_sl.rotateZ(side*60.*deg);
	      G4Transform3D rot_trans(rot_sl, G4ThreeVector(0., 0., 0.));
	      
	      G4Transform3D sl_back_trans
		= rot_trans*slice_back_transform[slice];

	      new G4PVPlacement(sl_back_trans, name+sl_name,
				slice_back_log[6*slice+side],
				det_env, false, 0);

	      sprintf(sl_name, "/%s_cap_slice_frontA_%d_%d",
		      name.c_str(), slice, side);
	      
	      G4Transform3D sl_frontA_trans
		= rot_trans*slice_frontA_transform[slice];

	      new G4PVPlacement(sl_frontA_trans, name+sl_name,
				slice_frontA_log[6*slice+side],
				det_env, false, 0);

	      if(slice_frontB_log[6*slice+side])
		{
		  sprintf(sl_name, "/%s_cap_slice_frontB_%d_%d",
			  name.c_str(), slice, side);
		  
		  G4Transform3D sl_frontB_trans
		    = rot_trans*slice_frontB_transform[slice];

		  new G4PVPlacement(sl_frontB_trans, name+sl_name,
				    slice_frontB_log[6*slice+side],
				    det_env, false, 0);

		}
	    }
	}
    }
  
  // install sensitive detector if HistoManager is set
  if(histoManager)
    {
      G4SDManager* SDman = G4SDManager::GetSDMpointer();
      HPGeSD = new MiniBallSensitiveDetector(histoManager, name);
      SDman->AddNewDetector(HPGeSD);
      HPGe_log->SetSensitiveDetector(HPGeSD);
    }
  
  return(mother);
}

// the segment number is stored as segment_id running from 0 to 5
// counterclockwise, for this is easier to handle.
// GetSegmentCode() gives back the number of the segment in the detector
      
G4int MiniBallSingleDetector::GetSegmentCode(G4int segment_nb)
{
  static G4int segmentCode[] = {3, 1, 2, 5, 6, 4};  

  if(segment_nb >= 0 && segment_nb < 6)
    return(segmentCode[segment_nb]);
  else
    return(-1);
}  

// returns the detectorID for the HistoManager

G4int MiniBallSingleDetector::GetDetectorID()
{
  if(HPGeSD)
    return(HPGeSD->GetDetectorID());
  else
    return(-1);
}

G4double MiniBallSingleDetector::GetDetectorEnvelopeLength()
{
  return(CAPSULE_HEIGHT + GAP_TO_CAPSULE + CAPSULE_FRONT_THICKNESS
	 + 2*SMALL_DIST);
}

G4double MiniBallSingleDetector::GetDetectorEnvelopeMidPosition()
{
  return(CAPSULE_HEIGHT/2+GAP_TO_CAPSULE+CAPSULE_FRONT_THICKNESS+SMALL_DIST);
}

G4double MiniBallSingleDetector::MoveToDetectorFront()
{
  return(CAPSULE_HEIGHT/2);
}

void MiniBallSingleDetector::SetDetectorVisAttributes(
  const G4VisAttributes* pvis_att)
{
  det_vis_att = pvis_att;
}

void MiniBallSingleDetector::SetDetectorVisAttributes(
  const G4VisAttributes& vis_att)
{
  det_vis_att = &vis_att;
}

void MiniBallSingleDetector::SetCapsuleVisAttributes(
  const G4VisAttributes* pvis_att)
{
  capsule_vis_att = pvis_att;
}

void MiniBallSingleDetector::SetCapsuleVisAttributes(
  const G4VisAttributes& vis_att)
{
  capsule_vis_att = &vis_att;
}
  

void MiniBallSingleDetector::CutCapsuleVisualisation(G4int type)
{
  G4cout << name << ": using cut capsule geometry "
	 << "- use this option only for visualisation!" << G4endl;
  capsule_cut=true;
  capsule_cut_type=type;
}

void MiniBallSingleDetector::CutDetectorVisualisation(G4int type)
{
  G4cout << name << ": using cut detector geometry "
	 << "- use this option only for visualisation!" << G4endl;
  det_cut=true;
  det_cut_type=type;
}

void MiniBallSingleDetector::SetEfficiencyFactor(G4double p_eff_fact)
{
  eff_fact = p_eff_fact;
}

void MiniBallSingleDetector::SetRelativeEfficiency(G4double percent)
{
  eff_fact = CalculateEfficiencyFactor(percent);
}

G4double MiniBallSingleDetector::CalculateEfficiencyFactor(G4double rel_eff)
{
  if(rel_eff>65.5 || rel_eff<16.)
    {
      G4cerr << "relative efficiency is out of bound, leaving unchanged"
	     << G4endl;
      return(1.);
    }

  G4double rel_eff2=rel_eff*rel_eff;
  G4double rel_eff3=rel_eff2*rel_eff;
  G4double rel_eff4=rel_eff3*rel_eff;

  return(EFF_REL_QUAD*rel_eff4+EFF_REL_THIRD*rel_eff3
    +EFF_REL_SECOND*rel_eff2+EFF_REL_LIN*rel_eff+EFF_REL_CONST);
}

void MiniBallSingleDetector::SetDeadlayer(G4double thickness, G4int type)
{
  deadlayer_thickness = thickness;
  deadlayer_type = type;
}
