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

#ifndef MiniBallSingleDetector_H
#define MiniBallSingleDetector_H 1

#include "G4VUserDetectorConstruction.hh"

// Include-files for detector construction

#include "globals.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"

#include "MBglobals.hh"
#include "MiniBallWorld.hh"
#include "MiniBallMaterial.hh"
#include "MiniBallGeom.hh"

// Inculde-files for detector readout
#include "G4SDManager.hh"
#include "MiniBallSensitiveDetector.hh"
#include "MiniBallHistoManager.hh"

#include "G4ios.hh"

// include-files for visualisation
#ifdef G4VIS_USE
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#endif

///Single HPGe crystal
class MiniBallSingleDetector : public G4VUserDetectorConstruction
{
public:
  /** \brief Constructor for a encapsulated MiniBall HPGe detector
   * 
   *  The position of the detector can be specified in various 
   *  ways, e.g. by polar coordinates rho, phi, theta or by a giving
   *  G4Transform3D. The parameter spin rotates the detector
   *  around its own axis. The distance rho from the origin is
   *  taken with respect to the front of the crystal.
   *  To position the detector with respect to the front
   *  of the capsule add the value returned by the static method
   *  MiniBallSingleDetector::MoveToDetectorFront().
   *  
   *  \param name name of the HPGe detector
   *  \param rho radius in polar coordinates
   *  \param theta theta angle in polar coordinates
   *  \param phi phi anlge in polar coordinates
   *  \param spin anlge for the rotation of the detector
   *         around its own axis
   *  \param transform transformation to position the triple-cluster
   *  \param mother mother volume, if this is not specified a mother
   *         volume is created (e.g. a MiniBallWorld object)
   *  \param histoManager a MiniBallHistoManager object for the readout
   *  \param capsule flag determining if the crystal capsule should be
   *         simulated
   *  \see MiniBallHistoManager
   *  \see MiniBallWorld
   */
  MiniBallSingleDetector(G4bool capsule=false);

  MiniBallSingleDetector(G4String name,
			 G4VPhysicalVolume* mother, G4bool capsule=false); 

  MiniBallSingleDetector(G4String name,
			 G4VPhysicalVolume* mother, G4double rho, 
			 G4double theta, G4double phi, G4bool capsule=false);

  MiniBallSingleDetector(G4String name,
			 G4VPhysicalVolume* mother, G4double rho, 
			 G4double theta, G4double phi, G4double spin, 
			 G4bool capsule=false);

  MiniBallSingleDetector(G4String name,
			 G4VPhysicalVolume* mother, G4Transform3D* transform, 
			 G4bool capsule=false); 

  MiniBallSingleDetector(G4String name,
			 G4VPhysicalVolume* mother, G4Transform3D* transform,
			 G4double spin, G4bool capsule=false); 

  MiniBallSingleDetector(G4double rho,
			 G4double theta, G4double phi, G4bool capsule=false);

  MiniBallSingleDetector(G4double rho, 
			 G4double theta, G4double phi, G4double spin, 
			 G4bool capsule=false);

  MiniBallSingleDetector(G4Transform3D* transform, 
			 G4bool capsule=false);

  MiniBallSingleDetector(G4Transform3D* transform,
			 G4double spin, G4bool capsule=false);

  MiniBallSingleDetector(MiniBallHistoManager* histoManager, 
			 G4bool capsule=false);

  MiniBallSingleDetector(MiniBallHistoManager* histoManager,
			 G4String name, G4VPhysicalVolume* mother, 
			 G4bool capsule=false); 

  MiniBallSingleDetector(MiniBallHistoManager* histoManager,
			 G4String name,G4VPhysicalVolume* mother, 
			 G4double rho, G4double theta, G4double phi, 
			 G4bool capsule=false);

  MiniBallSingleDetector(MiniBallHistoManager* histoManager,
			 G4String name, G4VPhysicalVolume* mother, 
			 G4double rho, G4double theta, G4double phi, 
			 G4double spin, G4bool capsule=false);

  MiniBallSingleDetector(MiniBallHistoManager* histoManager,
			 G4String name, G4VPhysicalVolume* mother,
			 G4Transform3D* transform, G4bool capsule=false); 

  MiniBallSingleDetector(MiniBallHistoManager* histoManager,
			 G4String name, G4VPhysicalVolume* mother, 
			 G4Transform3D* transform, G4double spin, 
			 G4bool capsule=false); 

  MiniBallSingleDetector(MiniBallHistoManager* histoManager,
			 G4double rho, G4double theta, G4double phi, 
			 G4bool capsule=false);

  MiniBallSingleDetector(MiniBallHistoManager* histoManager, 
			 G4double rho, G4double theta,G4double phi, 
			 G4double spin, G4bool capsule=false);
 
  MiniBallSingleDetector(MiniBallHistoManager* histoManager,G4Transform3D* transform, 
			 G4bool capsule=false);
  
  MiniBallSingleDetector(MiniBallHistoManager* histoManager,
			 G4Transform3D* transform,
			 G4double spin, G4bool capsule=false);

  ~MiniBallSingleDetector();

public:
  /** \brief construction of physical detector volumes
   *  
   *  To add a MiniBall detector in your own DetectorConstruction class
   *  create a MiniBallTripleCluster object and invoke the Construct()
   *  method.
   *
   *  \return mother volume
   */
  G4VPhysicalVolume* Construct();

  /** \brief returns the ID of the detector
   *
   *  The returned ID is needed to address the detector in the readout
   *
   *  \return ID of detector
   *  \see MiniBallHistoManager
   */
  G4int GetDetectorID();

  /** \brief returns the segment code of a segment specified by its ID
   *
   *  The segments of the detector are labeled with numbers on the cryostat
   *  in the following called segment code.
   *  In the simulation a segment ID is used running from 0 to 5 
   *  counterclockwise. This method returns the labeled segment code
   *  for a given segment specified by its ID.
   *
   *  \param seg_id segment ID
   *  \return segment code
   */
  static G4int GetSegmentCode(G4int seg_id);

  /** \brief return the length of the detector envelope volume
   *
   *  The detector is placed in an envelope volume. This method
   *  returns the length of this envelope volume. 
   *
   *  \return length of detector envelope volume
   */
  static G4double GetDetectorEnvelopeLength();

  /** \brief return the mid position of the detector envelope volume
   *
   *  The detector is placed in an envelope volume. This method
   *  returns the mid position of this envelope volume. 
   *
   *  \return mid position of detector envelope volume
   */
  static G4double GetDetectorEnvelopeMidPosition();

  /** \brief returns the negative distance from the crystal front
   *         to the front of the capsule 
   *
   *  The detector is positioned with respect to the front of the crystal
   *  To position the detector with respect to the front of the capsule
   *  add the returned value to your radius rho.
   *
   *  \return negative distance between from the crystal front to
   *          the front of the capsule
   */
  static G4double MoveToDetectorFront();

  /** \brief sets the relative efficiency of the detector
   *
   *  When setting the realtive efficiency of the detector the dimensions of
   *  the crystal are scaled correspondingly. 
   *  This is only a rough approximation.
   *  See also the efficiency example HPGeEff for more detail.
   *
   *  \param percent relative efficiency of the detector in percent
   *  \see HPGeEff
   */
  void SetRelativeEfficiency(G4double percent);

  /** \brief sets the scaling factor for the efficiency correction
   *         of the detector
   *
   *  When setting the scaling factor for the efficiency correction 
   *  of the detector the dimensions of the crystal are
   *  scaled correspondingly. This is only a rough approximation.
   *  See also the efficiency example HPGeEff for more detail.
   *
   *  \param eff_fact efficiency factor (scaling factor)
   *  \see HPGeEff
   */
  void SetEfficiencyFactor(G4double eff_fact);

  /** \brief approximates the efficiency factor for a given
   *         relative efficiency
   *
   *  When setting the realtive efficiency of the detector the dimensions of
   *  the crystal are scaled correspondingly. This method returns the
   *  scaling factor for a given relative efficiency.
   *  See also the efficiency example HPGeEff for more detail.
   *
   *  \param rel_eff relative efficiency of the detector in percent
   *  \return efficiency factor (scaling factor)
   *  \see HPGeEff
   */
  static G4double CalculateEfficiencyFactor(G4double rel_eff);

  /** \brief sets the thickness of the deadlayer
   *
   *  This method sets the thickness of the deadlayer and the type of
   *  deadlayer implementation (possible values are: DEADLAYER_FRONT - only
   *  the front of the detector has a deadlayer (default), 
   *  DEADLAYER_OUTER_CONTACT - the whole outer contact has a deadlayer
   *  (more realistic, but slower)
   *
   *  \param thickness thickness of the deadlayer
   *  \param type type of deadlayer implementation
   */
  void SetDeadlayer(G4double thickness, G4int type=DEFAULT_DEADLAYER_TYPE);
  
  /** \brief sets the visualisation attributes for the detector
   *
   *  \param pvis_att pointer to a G4VisAttributes object
   */
  void SetDetectorVisAttributes(const G4VisAttributes* pvis_att);

  /** \brief sets the visualisation attributes for the detector
   *
   *  \param vis_att G4VisAttributes object
   */
  void SetDetectorVisAttributes(const G4VisAttributes& vis_att);

  /** \brief sets the visualisation attributes for the detector capsule
   *
   *  \param pvis_att pointer to a G4VisAttributes object
   */  
  void SetCapsuleVisAttributes(const G4VisAttributes* pvis_att);

  /** \brief sets the visualisation attributes for the detector capsule
   *
   *  \param vis_att G4VisAttributes object
   */  
  void SetCapsuleVisAttributes(const G4VisAttributes& vis_att);

  /** \brief visualize a cut detector capsule
   *
   *  The geometry of the detector capsule is simulated cut open 
   *  to have a glance inside the structure of the simulated capsule.
   * 
   *  \warning Only use this method for visualisation means!
   *  \param type (0 for normal cut, 1 for cutting the other side)
   */
  void CutCapsuleVisualisation(G4int type=0);

  /** \brief visualize a cut detector
   *
   *  The geometry of the detector is simulated cut open 
   *  to have a glance inside the structure of the simulated detector.
   * 
   *  \warning Only use this method for visualisation means!
   *  \param type (0 for normal cut, 1 for cutting the other side)
   */
  void CutDetectorVisualisation(G4int type=0);

private:
  void GetTransformation();

private:
  G4String name;

  G4VPhysicalVolume* HPGe;
  G4VPhysicalVolume* deadlayer;
  G4VPhysicalVolume* Al_capsule; 
  G4VPhysicalVolume* det_env; 
  G4VPhysicalVolume* mother;
  MiniBallSensitiveDetector* HPGeSD;

  G4double rho;
  G4double phi;
  G4double theta;
  G4double spin;
  G4Transform3D* HPGe_transformation;

  G4bool capsule;

  MiniBallHistoManager* histoManager;

  G4double eff_fact;

  G4double deadlayer_thickness;
  G4int deadlayer_type;

  const G4VisAttributes* det_vis_att;
  const G4VisAttributes* capsule_vis_att;

  G4bool capsule_cut;
  G4bool det_cut;
  G4int capsule_cut_type;
  G4int det_cut_type;
};

#endif









