//////////////////////////////////////////////////////////////////////////////
/**  \file   MiniBallHistoManager.hh
 *   \brief  handling of the histogramming
 * 
 *   \author Hans Boie
 */
//////////////////////////////////////////////////////////////////////////////

#ifndef MiniBallHistoManager_h
#define MiniBallHistoManager_h 

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

#include "MiniBallHit.hh"
#include "MiniBallHistoGenerator.hh"
#include "MBglobals.hh"

///auxiliary class
class MiniBallDetectorHisto
{
public:
  MiniBallDetectorHisto(G4String);

  ~MiniBallDetectorHisto();
  
  G4String DetectorName;
  G4int DetectorHitsCollectionID;
  
  G4int cryostatID;   
  
  MiniBallHistogram* CoreHisto;
  G4int CoreEnNBin;
  G4double CoreEnLow, CoreEnHigh;
  G4double CoreEn, CoreEnAccu;
  G4double time;

  G4double MIEnergy;
  G4int MIHitSegment;
  G4ThreeVector MIPositionInDet;
  G4ThreeVector TRInteractionPosition;	

  MiniBallHistogram* SegmentHisto[NB_OF_SEGMENTS];
  G4int SegmentEnNBin[NB_OF_SEGMENTS];
  G4double SegmentEnLow[NB_OF_SEGMENTS], SegmentEnHigh[NB_OF_SEGMENTS];
  G4double SegmentEn[NB_OF_SEGMENTS], SegmentEnAccu[NB_OF_SEGMENTS];

  G4int MBTupleCoreEnId;
  G4int MBTupleSegmentEnId[NB_OF_SEGMENTS];

  G4double CoreElectronicSigma;
  G4double CoreStatisticEpsilon;
  G4double SegmentElectronicSigma[NB_OF_SEGMENTS];
  G4double SegmentStatisticEpsilon[NB_OF_SEGMENTS];

  G4double CoreCutoffEn;
  G4double CoreCutoffWidth;
  G4double CoreSuppressionEn;
  G4double CoreSuppressionSlope;
  G4double SegmentCutoffEn[NB_OF_SEGMENTS];
  G4double SegmentCutoffWidth[NB_OF_SEGMENTS];
  G4double SegmentSuppressionEn[NB_OF_SEGMENTS];
  G4double SegmentSuppressionSlope[NB_OF_SEGMENTS];
}; 

/// auxiliary class
class MiniBallCryostatHisto
{
public:
  MiniBallCryostatHisto(G4String);
  ~MiniBallCryostatHisto();

  G4String CryostatName;
   
  MiniBallHistogram* CryostatHisto;
  G4int CryostatEnNBin;
  G4double CryostatEnLow, CryostatEnHigh;

  G4int MBTupleCryostatEnId;
  std::vector<size_t> detectorIDs;
}; 

/// auxiliary class
class MiniBallArray
{
public:
  MiniBallArray() {;};
  ~MiniBallArray() {;};

  G4int MBTupleArrayID;
  std::vector<G4int> segments;

  G4int MBTupleHitSegmentID;
  G4int MBTupleRadiusInDetID;
  G4int MBTuplePhiInDetID;
  G4int MBTupleZInDetID;
  G4int MBTupleXInDetID;
  G4int MBTupleYInDetID;
};

class MiniBallHistoManager 
{
public:
  MiniBallHistoManager(MiniBallHistoGenerator* histo_generator);
  ~MiniBallHistoManager();

public:
  
//////////////////////////////////////////////////////////////////////////////
///\name Creating histogramms
//////////////////////////////////////////////////////////////////////////////
//@{

  /// Generate a histogram for a given segment in a certain detector.
  /** The detector can be spezified by its cryostat and its number in 
      that cryostat (det_nb), or with its absolut number (det_id). 
      The core can be referred to as the segment with the id CORE. 
      The name of the histogramm can be specified.*/
  /// det_nb is the number of the detector within a given cryostat,
  /// it can have the values DET_A, DET_B, DET_C 

  void BookHisto(G4String name, G4int cryo_id, G4int det_nb, G4int seg_id);
  /// \overload
  void BookHisto(G4int cryo_id, G4int det_nb, G4int seg_id);

  /// \overload  
  /// det_id is the absolute detector id
  void BookHisto(G4String name, G4int det_id, G4int seg_id);
  /// \overload
  void BookHisto(G4int det_id, G4int seg_id);
  
  /// \overload
  void BookHisto(G4String name,G4int seg_id);
  /// \overload
  void BookHisto(G4int seg_id);
  
////////////////////////////////////////////////////////////////////////////  

  /// Generate a histograms for the addback-energy in a given cryostat.
  void BookAddbackCryostat(G4String name, G4int cryo_id);
  /// \overload
  void BookAddbackCryostat(G4int cryo_id);

//////////////////////////////////////////////////////////////////////////////
  
  /// Generate all histograms in a given detector.
  /// det_nb is the number of the detector within a given cryostat,
  /// it can have the values DET_A, DET_B, DET_C 
  void BookAllInDetector(G4String name, G4int cryo_nb, G4int det_id);
  /// \overload
  void BookAllInDetector(G4int cryo_id, G4int det_nb);

  /// \overload
  /// det_id is the absolute detector id
  void BookAllInDetector(G4String name, G4int det_id);
  /// \overload
  void BookAllInDetector(G4int det_id);

  /// Generate all histograms for all segments of all detectors in a cryostat
  // (including cores) and for the addback energy of the cryostat.
  void BookAllInCryostat(G4String name, G4int cryo_id);
  /// \overload
  void BookAllInCryostat(G4int cryo_id);
    
  // Generate histograms for all segments of all detectors in all cryostat
  // (including cores) and for the addback energys of the cryostats.
  void BookAll(G4String name);
  /// \overload
  void BookAll();

//@}
//////////////////////////////////////////////////////////////////////////////
///\name Creating Tuple-Fields
//////////////////////////////////////////////////////////////////////////////
//@{

  /// Add field for the energy in a given segment of a certain detector.
  /// det_nb is the number of the detector within a given cryostat,
  /// it can have the values DET_A, DET_B, DET_C 
  void AddFieldSegmentEn(G4String name, G4int cryo_id, G4int det_nb,
			 G4int seg_id);
  ///\overload  
  void AddFieldSegmentEn(G4int cryo_id, G4int det_nb, G4int seg_id);

  /// Add field for the energy in a given segment of a certain detector.
  /// det_id is the absolute detector id
  void AddFieldSegmentEn(G4String name, G4int det_id, G4int seg_id);
  ///\overload
  void AddFieldSegmentEn(G4int det_id, G4int seg_id);

  /// Add field for the energy in a given segment of the detector with
  /// detector id NULL
  void AddFieldSegmentEn(G4String name, G4int seg_id);
  ///\overload
  void AddFieldSegmentEn(G4int seg_id);

  // Add field for the addback energy of a given cryostat.
  void AddFieldAddbackCryostatEn(G4String name, G4int cryo_id);
  ///\overload
  void AddFieldAddbackCryostatEn(G4int cryo_id);

///////////////////////////////////////////////////////////////////////////////

  /// Add field for the energies of all segments in a certain detector.
  /// det_nb is the number of the detector within a given cryostat,
  /// it can have the values DET_A, DET_B, DET_C 
  void AddFieldAllInDetectorEn(G4String name,G4int cryo_id, G4int det_nb);
  ///\overload
  void AddFieldAllInDetectorEn(G4int cryo_id, G4int det_nb);

  ///\overload
  /// det_id is the absolute detector id
  void AddFieldAllInDetectorEn(G4String name, G4int det_id);
  ///\overload
  void AddFieldAllInDetectorEn(G4int det_id);
    
  /// Add field for the energies of all segments of all detectors of
  /// a given cryostat and the addback energy of the cryostat.
  void AddFieldAllInCryostatEn(G4String name, G4int cryo_id);
  ///\overload
  void AddFieldAllInCryostatEn(G4int cryo_id);

  /// Add field for the energies of all segments of all detectors of
  /// all cryostats and the addback energies of the cryostats.
  void AddFieldAllEn(G4String name);
  ///\overload
  void AddFieldAllEn();

/////////////////////////////////////////////////////////////////////////////

  ///Add field for the number of the hit segment, running from 0 to 5 clockwise
  void AddFieldSegmentNb(G4String name="SegmentNb");

  ///Add field for the real number (as written on the cryostat) of the
  ///hit segment.
  void AddFieldSegmentCode(G4String name="SegmentCode");

  /// Add field for the id of the hit detector.
  void AddFieldDetectorNb(G4String name="DetectorNb");

  /// Add field for the id of the hit cryostat. 
  void AddFieldCryostatNb(G4String name="CryostatNb");

  /// Add field for the number of the event in which this hit occurs
  void AddFieldEventNb(G4String name="EventNb");
 
  /// Add field for the trace,e.g. after which fundamental interaction 
  void AddFieldTraceNb(G4String name="TraceNb");

//////////////////////////////////////////////////////////////////////////////

  /// Add fields (r, phi and z) for the position in the hit detector.
  void AddFieldPositionInDet(G4String name= "Det");
   
  // Add fields r, phi and z separately.
  ///Add field r in detector
  void AddFieldRadiusInDet(G4String name="DetR");
  ///Add field phi in detector
  void AddFieldPhiInDet(G4String name="DetPhi");
  ///Add field z in detector
  void AddFieldZInDet(G4String name="DetZ");

  ///Add fields x
  void AddFieldXInDet(G4String name="DetX");
  ///Add fiel y 
  void AddFieldYInDet(G4String name="DetY");

/////////////////////////////////////////////////////////////////////////////  

  /// Add fields (world_x, world_y and world_z) for the position 
  /// of the hit in the world volume.
  void AddFieldPositionInWorld(G4String name="World");

  // Add fields world_x, world_y and world_z separately. 
  ///Add field worldX
  void AddFieldWorldX(G4String name="WorldX");
  ///Add field worldY
  void AddFieldWorldY(G4String name="WorldY");
  ///Add field WorldZ
  void AddFieldWorldZ(G4String name="WorldZ");
  ///Add field WorldTheta
  void AddFieldWorldTheta(G4String name="WorldTheta");
  ///Add field WorldPhi
  void AddFieldWorldPhi(G4String name="WorldPhi");

////////////////////////////////////////////////////////////////////////////  

  /// Add fields (r, phi and z) for the position of the first interation 
  ///in the hit detector.
  void AddFieldFIPositionInDet(G4String name="FIDet");

  // Add fields r, phi and z separately.
  ///Add field r of first interaction
  void AddFieldFIRadiusInDet(G4String name="FIDetR");
  ///Add field phi of first interaction
  void AddFieldFIPhiInDet(G4String name="FIDetPhi");
  ///Add field z of first interaction
  void AddFieldFIZInDet(G4String name="FIDetZ");

  ///Add fields x of first interaction
  void AddFieldFIXInDet(G4String name="FIDetX");
  ///Add field y of first interaction
  void AddFieldFIYInDet(G4String name="FIDetY");

///////////////////////////////////////////////////////////////////////////////

  /// Add fields (world_x, world_y and world_z) for the position 
  /// of the first interaction in the world volume.
  void AddFieldFIPositionInWorld(G4String name="FIWorld");

  //name Add fields world_x, world_y and world_z separately.  
  ///Add field FIWorldX 
  void AddFieldFIWorldX(G4String name="FIWorldX");
  ///Add field FIWorldY
  void AddFieldFIWorldY(G4String name="FIWorldY");
  ///Add field FIWorldZ
  void AddFieldFIWorldZ(G4String name="FIWorldZ");
  ///Add field FIWorldTheta
  void AddFieldFIWorldTheta(G4String name="FIWorldTheta");
  ///Add field FIWorldPhi
  void AddFieldFIWorldPhi(G4String name="FIWorldPhi");

//////////////////////////////////////////////////////////////////////////////

  /// Add user defined array. This allows you to generate an output-file
  /// similar to the output of the DAQ.
  void AddUserDefinedArray(G4String name,
			   G4int nb_of_entries, G4int* segment_nb);

  /// Array with all core energies
  void AddCoresArray(G4String name="Core");

  /// Adding array with hit segment number or position of main interaction
  void AddArrayHitSegment(G4String name="HitSeg");
  void AddArrayPositionInDet(G4String name="Det");
  void AddArrayRadiusInDet(G4String name="DetR");
  void AddArrayPhiInDet(G4String name="DetPhi");
  void AddArrayZInDet(G4String name="DetZ");
  void AddArrayXInDet(G4String name="DetX");
  void AddArrayYInDet(G4String name="DetY");
//@}
//////////////////////////////////////////////////////////////////////////////
///\name Internal organization
//////////////////////////////////////////////////////////////////////////////
//@{
  
  //The following methods are for internal use and should not be invoked
  //by the user directly.
  
  /// Add a detector to the list of registered detectors. The detector
  /// name can be specified. The detector id will be returned.
  G4int AddDetector(G4String name);
  ///\overload
  G4int AddDetector();
  
  /// Add a cryostat to the list of registered cryostat. 
  /** The cryostat name can be specified. The cryostat id will be returned.*/
  G4int AddCryostat(G4String name);
  ///\overload
  G4int AddCryostat();

  /// Connect a given detector to a certain cryostat. 
  /**The detector and 
     cryostat must be referred to via the detector or cryostat id returned
     by the methods AddMBDetector() or AddMBCryostat() respectively.*/
  void AddToCryostat(G4int cryo_id, G4int det_id);
  
//////////////////////////////////////////////////////////////////////////////

  /// Returns the number of registered detectors. 
  G4int GetNbOfDetectors();

  /// Returns the number of registered cryostats.
  G4int GetNbOfCryostats();

  // included by S.Klupp
  std::vector<MiniBallDetectorHisto*>* GetDetectors() {return &detectors; };
  std::vector<MiniBallCryostatHisto*>* GetCryostats() {return &cryostats; };

///////////////////////////////////////////////////////////////////////////////

  /// Sets the HitsCollection.
  void SetMBHitsCollectionID(G4int det_id, G4int HCID);
  
  /// Returns the HitsCollection ID.
  G4int GetMBHitsCollectionID(G4int det_id);

//////////////////////////////////////////////////////////////////////////////

  /// Set the binning and range of the histogram for a given segment 
  /// in a certain detector. The core can be referred to as segment CORE.
  void SetMBEnNBin(G4int det_id, G4int seg_id, G4int NBin);
  void SetMBEnLow(G4int det_id, G4int seg_id, G4double Low);
  void SetMBEnHigh(G4int det_id, G4int seg_id, G4double High);

  /// Set the binning and range of all booked histograms for a given detector. 
  void SetMBEnNBin(G4int det_id, G4int NBin);
  void SetMBEnLow(G4int det_id, G4double Low);
  void SetMBEnHigh(G4int det_id, G4double High);

  // Set the binning and range of the addback histogram for a given cryostat.
  void SetMBCryostatEnNBin(G4int cryo_id, G4int NBin);
  void SetMBCryostatEnLow(G4int cryo_id, G4double Low);
  void SetMBCryostatEnHigh(G4int cryo_id, G4double High);

  /// Set the binning and range for all booked histograms
  void SetMBEnNBin(G4int NBin);
  void SetMBEnLow(G4double Low);
  void SetMBEnHigh(G4double High);

//////////////////////////////////////////////////////////////////////////////

  /// Set energy unit for histograms and energy fields.
  /// Use energy units as defined in Geant4, like MeV, keV,...
  void SetEnergyUnit(G4double energy_unit);

//////////////////////////////////////////////////////////////////////////////

  /// Set sigma for gaussian line shape
  /// the electronic noise is added independently from the energy
  /// <s_E> = <el_sigma>
  /// the noise due to statistic fluctuation in the number of cahrge carriers
  /// goes with the square root of the detected energy
  /// <s_D> = sqrt(<stat_epsilon>*<energy>)
  /// the resulting resolution is given by
  /// <s_T> = sqrt(<s_E>^2 + <s_D>^2)
  void SetSigma(G4int det_id, G4int seg_id, G4double el_sigma=DEFAULT_EL_SIGMA,
		G4double stat_epsilon=DEFAULT_STAT_EPSILON);
  void SetAllSigma(G4int seg_id, G4double el_sigma=DEFAULT_EL_SIGMA,
		   G4double stat_epsilon=DEFAULT_STAT_EPSILON);
  void SetAllSigma(G4double el_sigma=DEFAULT_EL_SIGMA,
		   G4double stat_epsilon=DEFAULT_STAT_EPSILON);

  void SetThreshold(G4int det_id, G4int seg_id,
		    G4double cutoff_en=DEFAULT_CUTOFF_EN,
		    G4double cutoff_width=DEFAULT_CUTOFF_WIDTH,
		    G4double suppresion_en=DEFAULT_SUPPRESSION_EN,
		    G4double suppresion_slop=DEFAULT_SUPPRESSION_SLOPE);
  void SetAllThreshold(G4int seg_id, G4double cutoff_en=DEFAULT_CUTOFF_EN,
		       G4double cutoff_width=DEFAULT_CUTOFF_WIDTH,
		       G4double suppresion_en=DEFAULT_SUPPRESSION_EN,
		       G4double suppresion_slop=DEFAULT_SUPPRESSION_SLOPE);
  void SetAllThreshold(G4double cutoff_en=DEFAULT_CUTOFF_EN,
		       G4double cutoff_width=DEFAULT_CUTOFF_WIDTH,
		       G4double suppresion_en=DEFAULT_SUPPRESSION_EN,
		       G4double suppresion_slop=DEFAULT_SUPPRESSION_SLOPE);

//////////////////////////////////////////////////////////////////////////////

  // The following methods handle the creation of an ntuple for the
  // hits in the detector.

  /** SetMBTupleType set the type of tuple which should be used. You
      can either get one entry per Event (EVENT), so the accumulated energies
      of the event will be stored in the tuple and the detector number,
      segment number and position of the main interaction (MI) will be used,
      or get an entry for every interaction (INTERACTION), then only the 
      energie of this interaction, detector and segment numbers and position 
      will be stroed in the ntuple.*/
  void SetMBTupleType(G4int type) { MBTupleType=type; };

  /// Set name of ntuple.
  void SetMBTupleName(G4String tuple_name) { MBTupleName=tuple_name; };

  /// Get name of ntuple.
  G4String GetMBTupleName() { return MBTupleName; };

///////////////////////////////////////////////////////////////////////////////

  // Set or get the MiniBallTuple used. This method can be used to
  // add user-specified fields to the ntuple via the AddField() method
  // of the class MiniBallTuple. To do so the user can either generate
  // his own MiniBallTuple and register it for the MiniBall fields using
  // the method SetMiniBallTuple() or add his own fields to the
  // automatically generated ntuple returned by the method GetMiniBallTuple().

  // The method GetMiniBallTuple() must be invoked after generating
  // a least one MiniBall field by usind one of the AddMBField...() methods.
  MiniBallTuple* GetMiniBallTuple();
  

  /** The method SetMiniBallTuple() must be invoked before calling any
      of the AddMBField...() methods. If called afterwards it will not show
      any effect.*/
  void SetMiniBallTuple(MiniBallTuple* tup);

  /// Check if there were any data writen by the MiniBall Clusters. 
  /** (If the the MiniBallTuple is managed by the user, it can be checked in
    order not to write empty ntuples)*/
  G4bool WriteMBTuple() { return(WroteMBTuple); }

  /// Report that data has been written into the ntuple.
  void SetMBTupleWritten() { WroteMBTuple=true; };

//////////////////////////////////////////////////////////////////////////////

  /// Set position of the interaction.
  void SetMBInteraction(size_t det_id, size_t segment_id, G4double en,
			G4int event_id,G4int trace_id,
			G4ThreeVector positionInDet,
			G4ThreeVector positionInWorld,
			G4double time);

//////////////////////////////////////////////////////////////////////////////

  /// Initialize accumulation of energies
  void StartOfEvent();

  /// Fill histograms and ntuple if MBTupleType is EVENT.
  void EndOfEvent();

  /// Clear Interaction values
  void StartOfInteraction();

  /// Fill ntuple if MBTupleType is INTERACTION.
  void EndOfInteraction();
  //@}

  /// Get accumulated energy in a given segment of a given detector
  /// in a given cryostat
  /// det_nb is the relative detector number in the cryostat
  /// it can have the values DET_A, DET_B, DET_C 
  G4double GetSegmentEn(size_t cryo_id, size_t det_nb, size_t segment_id);

  /// Get accumulated energy in a given segment of a given detector
  /// det_id is the absolute detector id
  G4double GetSegmentEn(size_t det_id, size_t segment_id);

  /// Get accumulated energy of a given cryostat
  G4double GetAddbackEn(size_t cryo_id);

  //protected:
  /// Initialize accumulation of energy during an event.
  void InitMBAccumulation();

  /// Initialize values for new interaction
  void InitMBInteraction();

  /// Add gaussian lineshape
  void AddGaussian();

protected :

  /// Add gaussian lineshape to given energy
  G4double RandSigma(G4double eletronic_sigma, G4double statistic_sigma,
		     G4double energy);

  /// Cut off energy spectra
  void SetCutoff();

  /// Calculate energy cutoff shape
  G4bool RandCutoff(G4double en, G4double cutoff_en, G4double cutoff_width,
		    G4double suppession_en, G4double suppression_slope);

  /// Write accumulated energies in histogramms.
  void FillMBHisto();

  /// Fill the ntuple at the end of an event.
  void FillMBTuple(); 

  /// Fill the ntuple at the end of an interaction.
  void FillMBInteractionTuple(); 

  MiniBallHistoGenerator* HistoGenerator;

  std::vector<MiniBallDetectorHisto*> detectors;
  std::vector<MiniBallCryostatHisto*> cryostats;

  std::vector<MiniBallArray*> arrays;

  MiniBallTuple* MBTuple;
  
  G4String MBTupleName;
  G4int MBTupleType;
  G4int MBTupleManagement;

  G4int MBTupleSegmentNbId, MBTupleSegmentCodeId;
  G4int MBTupleDetectorNbId, MBTupleCryostatNbId;
  G4int MBTupleEventNbId;
  G4int MBTupleTraceNbId;
  G4int MBTupleRadiusId, MBTuplePhiId, MBTupleZId;
  G4int MBTupleXId, MBTupleYId;

  G4int MBTupleWorldXId, MBTupleWorldYId, MBTupleWorldZId;
  G4int MBTupleWorldThetaId, MBTupleWorldPhiId;

  G4int MBTupleFIRadiusId, MBTupleFIPhiId, MBTupleFIZId;
  G4int MBTupleFIXId, MBTupleFIYId;

  G4int MBTupleFIWorldXId, MBTupleFIWorldYId, MBTupleFIWorldZId;
  G4int MBTupleFIWorldThetaId, MBTupleFIWorldPhiId;

  G4bool AnyMBTupleCryostatEn;
  G4bool AnyMBTupleSegmentEn;
  G4bool AnyMBTupleCoreEn;

  G4bool WroteMBTuple;
  G4bool Interaction;
  G4bool InteractionInEvent;

  G4int InteractionDetectorNb, InteractionSegmentNb, InteractionCryostatNb;
  G4double InteractionEn;
  G4int InteractionEventNb;
  G4int InteractionTraceNb;

  G4ThreeVector InteractionPositionInDet;
  G4ThreeVector InteractionPositionInWorld;

  G4int MainInteractionDetectorNb, MainInteractionSegmentNb;
  G4int MainInteractionCryostatNb;
  G4double MainInteractionEn;

  G4ThreeVector MainInteractionPositionInDet;
  G4ThreeVector MainInteractionPositionInWorld;

  G4ThreeVector FirstInteractionPositionInDet;
  G4ThreeVector FirstInteractionPositionInWorld;

  G4double EnergyUnit;

  G4bool WriteGaussian;
  G4bool ElectronicThreshold;
};

#endif










