/************************************************************************
 * \file MBglobals.hh
 *
 * \brief global definitions used in the MiniBall package
 *
 * All definitions can easily be varied without changing this file
 * by using the "-D" option of gcc.
 * You can use the environment variable MB_DEFINITIONS, e.g. put the
 * following line into your makefile: 
 * 
 *     export MB_DEFINITIONS := -DGAP_BETWEEN_DETECTORS=1.5 
 *
 * into your makefile to set the definition of GAP_BETWEEN_DETECTORS to 
 * 1.5mm
 *
 * \author Tanja Striepling (tanja@ikp.uni-koeln.de)
 * \author hans.boie@mpi-hd.mpg.de
 * \author v.bildstein@mpi-hd.mpg.de
 *
 ************************************************************************/

#ifndef MBGLOBALS_H
#define MBGLOBALS_H 

#include "globals.hh"
#include "G4Material.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

/************************************************************************
 * \name Crystall-dimensions:
 *
 ************************************************************************/
//@{

///height of the detector-crystall  
#ifndef DETECTOR_HEIGHT
#define DETECTOR_HEIGHT (77.71*mm)
#endif
///diameter: corner to corner at the hexagonal front
#ifndef DETECTOR_FRONT
#define DETECTOR_FRONT (67.89*mm)
#endif
///diameter: edge to edge at the hexagonal front
#ifndef SW_DETECTOR_FRONT
#define SW_DETECTOR_FRONT (DETECTOR_FRONT*sqrt(3.)/2.)
#endif
///diameter: at the back of the detector
#ifndef DETECTOR_BACK
#define DETECTOR_BACK (70.00*mm)
#endif
///front radius of the tub, which helps to compose the detectors shape
#ifndef DETECTOR_MIN_RADIUS
#define DETECTOR_MIN_RADIUS (SW_DETECTOR_FRONT/2.)
#endif
///back radius of the tub, which helps to compose the detectors shape
#ifndef DETECTOR_MAX_RADIUS
#define DETECTOR_MAX_RADIUS (DETECTOR_BACK/2.)
#endif
///angle between the front plate and the sides
#ifndef DETECTOR_ANGLE
#define DETECTOR_ANGLE (4.125*deg)
#endif
///thickness of deadlayer of detector
#ifndef DETECTOR_DEFAULT_DEADLAYER_THICKNESS
#define DETECTOR_DEFAULT_DEADLAYER_THICKNESS (0.0005*mm)
#endif
///posible implementations of the deadlayer of the crystal
#define DEADLAYER_FRONT 0
#define DEADLAYER_OUTER_CONTACT 1

#ifndef DEFAULT_DEADLAYER_TYPE
#define DEFAULT_DEADLAYER_TYPE DEADLAYER_FRONT
#endif

///gap between the crystal and its capsule
#ifndef GAP_TO_CAPSULE
#define GAP_TO_CAPSULE (0.7*mm)
#endif
///thickness of the aluminium capsule
#ifndef CAPSULE_THICKNESS
#define CAPSULE_THICKNESS (0.7*mm)
#endif
#ifndef CAPSULE_FRONT_THICKNESS
#define CAPSULE_FRONT_THICKNESS (1.*mm)
#endif
///length of the capsule
#ifndef CAPSULE_HEIGHT
#define CAPSULE_HEIGHT (100.*mm)
#endif
///position of capsule back
#ifndef CAPSULE_BACK
#define CAPSULE_BACK (94.*mm)
#endif
///radius of contact hole
#ifndef CAPSULE_CONTACT_HOLE_RAD
#define CAPSULE_CONTACT_HOLE_RAD (8.75*mm)
#endif
///depth of contact hole
#ifndef CAPSULE_CONTACT_HOLE_DEPTH
#define CAPSULE_CONTACT_HOLE_DEPTH (4.5*mm)
#endif

///number of segments
#ifndef NB_OF_SEGMENTS
#define NB_OF_SEGMENTS 6
#endif

///the solid of the cap is approximated as build up out of
///flat slices (at least 7, must be odd)
#ifndef NB_OF_SLICES
#define NB_OF_SLICES 5
#endif

///very small distance
#ifndef SMALL_DIST
#define SMALL_DIST (0.01*mm)
#endif
//@}

/************************************************************************
 * \name Cryostat-dimensions
 *
 ************************************************************************/
//@{

/// 
#ifndef GAP_TO_CRYOSTAT
#define GAP_TO_CRYOSTAT (2.0*mm)
#endif
/// thickness of the wall of the cyostat
#ifndef CRYOSTAT_THICKNESS
#define CRYOSTAT_THICKNESS (1.5*mm)
#endif
/// gap between the capsules in the cryostat
#ifndef GAP_BETWEEN_DETECTORS
#define GAP_BETWEEN_DETECTORS (0.5*mm)
#endif
/// 
#ifndef CRYOSTAT_FRONT
#define CRYOSTAT_FRONT ( (SW_DETECTOR_FRONT \
			  +2*(GAP_TO_CAPSULE+CAPSULE_THICKNESS) \
			  *(1/cos(DETECTOR_ANGLE)-tan(DETECTOR_ANGLE))) \
			 *cos(asin(2/sqrt(3.)*sin(DETECTOR_ANGLE))) \
                        +GAP_BETWEEN_DETECTORS+GAP_TO_CRYOSTAT \
                        +CRYOSTAT_THICKNESS)
#endif
///
#ifndef CAP_RADIUS
#define CAP_RADIUS (114.42*mm)
#endif
///
#ifndef CAP_HEIGHT
#define CAP_HEIGHT (4.*mm)
#endif
///
#ifndef CAP_DELTA_THETA
#define CAP_DELTA_THETA (15.19*deg)
#endif
///
const G4double CRYOSTAT_SIDE_HEIGHT1=91.5*mm;
///
const G4double CRYOSTAT_SIDE_HEIGHT2=35.*mm;
///
const G4double CRYOSTAT_SIDE_THICKNESS2=8.*mm;
/// 
const G4double CRYOSTAT_NOZZLE_INNER_RAD[]={31.*mm, 31.*mm, 29.*mm, 29.*mm,
                 29.*mm, 29.*mm, 28.5*mm, 28.5*mm, 29.5*mm, 29.5*mm};
///
const G4double CRYOSTAT_NOZZLE_OUTER_RAD[]={35.*mm, 35.*mm, 35.*mm, 35.*mm, 
		 32.5*mm, 32.5*mm, 32.5*mm, 32.5*mm, 32.5*mm, 32.5*mm};
///
const G4double CRYOSTAT_NOZZLE_Z_PLANE[]={0.*mm, -9.*mm, -11.*mm, -12.*mm,
		 -12.*mm, -18.5*mm, -19.*mm, -24.5*mm, -41.*mm};

const G4int CRYOSTAT_NOZZLE_NB_OF_Z_PLANES=9;
/// Length of the world volume
#ifndef WORLD_LENGTH
#define WORLD_LENGTH (3*m)
#endif
/**< actually the world ist a box, which is twice the WORLD_LENGTH
   big in all three directions*/

/// cluster types
#define TRIPLE_CLUSTER 3
//@}

/************************************************************************
 * \name Histograming
 *
 ************************************************************************/
//@{

///detector in cryostat
#define DET_A 0
///
#define DET_B 1
///
#define DET_C 2
///
#define DET_D 3

///segments
#define SEG_0 0
///
#define SEG_1 1
///
#define SEG_2 2
///
#define SEG_3 3
///
#define SEG_4 4
///
#define SEG_5 5

///segment number of the core
#define CORE NB_OF_SEGMENTS+1
///
#define ALL_SEG CORE+1
///
#define EMPTY_CHN ALL_SEG+1

///id for cryostats
#define CRYOSTAT -2

///flag which indicates if there should be an tuple created for this field
#define FIELD_NOT_SET -1

///flag indicating wether the filling is done automatically or not
#define EXTERN_MANAGEMENT 0
///
#define INTERN_MANAGEMENT 1

///flag for the type of MiniBall tuple
///all interactions in one event are summed up to one row in the tuple
#define EVENT 0
///each interaction is one row in the tuple
#define INTERACTION 1

///default binning for the histogramms
#ifndef DEFAULT_ENNBIN
#define DEFAULT_ENNBIN 16000
#endif
///the lower end of the range in the histogramms
#ifndef DEFAULT_ENLOW
#define DEFAULT_ENLOW 0.
#endif
///the higher end of the range in the histogramms
#ifndef DEFAULT_ENHIGH
#define DEFAULT_ENHIGH (3.2)
#endif
///energy unit of histogramming
#ifndef DEFAULT_ENERGY_UNIT
#define DEFAULT_ENERGY_UNIT MeV
#endif

/// way of collecting data
#ifndef MB_TUPLE_TYPE
#define MB_TUPLE_TYPE EVENT
#endif

///default tuple name
#ifndef DEFAULT_TUPLE_NAME
#define DEFAULT_TUPLE_NAME "MB"
#endif

///name of hits collection
#ifndef MB_HITS_COLLECTION_NAME
#define MB_HITS_COLLECTION_NAME "MBHitsCollection"
#endif

///detector resolution
#define DEFAULT_EL_SIGMA (1.3*keV)
#define DEFAULT_STAT_EPSILON (0.*keV)

///energy threshold
#define DEFAULT_CUTOFF_EN (30.*keV)
#define DEFAULT_CUTOFF_WIDTH (5.*keV)
#define DEFAULT_SUPPRESSION_EN 0.
#define DEFAULT_SUPPRESSION_SLOPE -1.
//@}

/************************************************************************
 * \name Relative efficiency relation
 *
 ************************************************************************/
//@{

#define EFF_REL_QUAD -1.69e-08
#define EFF_REL_THIRD 3.50e-06
#define EFF_REL_SECOND -0.000304
#define EFF_REL_LIN 0.0176
#define EFF_REL_CONST 0.482

//@}

#endif
