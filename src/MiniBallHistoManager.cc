/************************************************************************
 * 
 *  handle histogramming
 * 
 ************************************************************************/

#include "MiniBallHistoManager.hh"

#include "MBglobals.hh"
#include "MiniBallSingleDetector.hh"

#include "G4UnitsTable.hh"
#include <string>

MiniBallDetectorHisto::MiniBallDetectorHisto(G4String name)
  :DetectorName(name),DetectorHitsCollectionID(-1),cryostatID(-1),
   CoreHisto(NULL),CoreEnNBin(DEFAULT_ENNBIN),CoreEnLow(DEFAULT_ENLOW),
   CoreEnHigh(DEFAULT_ENHIGH),CoreEn(0),
   CoreEnAccu(0),MBTupleCoreEnId(FIELD_NOT_SET)
{
  for(int i=0; i<NB_OF_SEGMENTS;i++)
    {
      SegmentHisto[i]=NULL;
      SegmentEnNBin[i]=DEFAULT_ENNBIN;
      SegmentEnLow[i]=DEFAULT_ENLOW;
      SegmentEnHigh[i]=DEFAULT_ENHIGH;
      SegmentEn[i]=0;
      SegmentEnAccu[i]=0;
      MBTupleSegmentEnId[i]=FIELD_NOT_SET;
    }
}

MiniBallCryostatHisto::MiniBallCryostatHisto(G4String name)
  :CryostatName(name),CryostatHisto(NULL),CryostatEnNBin(DEFAULT_ENNBIN),
   CryostatEnLow(DEFAULT_ENLOW), CryostatEnHigh(DEFAULT_ENHIGH),
   MBTupleCryostatEnId(FIELD_NOT_SET)
{;}

MiniBallHistoManager::MiniBallHistoManager(
                           MiniBallHistoGenerator* pHistoGenerator)
  : HistoGenerator(pHistoGenerator), detectors(0), cryostats(0), MBTuple(NULL),
    MBTupleName(DEFAULT_TUPLE_NAME), MBTupleType(MB_TUPLE_TYPE),
    MBTupleManagement(INTERN_MANAGEMENT), AnyMBTupleCryostatEn(false),
    AnyMBTupleSegmentEn(false), AnyMBTupleCoreEn(false),
    EnergyUnit(DEFAULT_ENERGY_UNIT), WriteGaussian(false),
    ElectronicThreshold(false)
{
  MBTupleSegmentNbId=FIELD_NOT_SET;
  MBTupleSegmentCodeId=FIELD_NOT_SET;
  MBTupleDetectorNbId=FIELD_NOT_SET;
  MBTupleCryostatNbId=FIELD_NOT_SET;
  MBTupleEventNbId=FIELD_NOT_SET;
  MBTupleTraceNbId=FIELD_NOT_SET;
  MBTupleRadiusId=FIELD_NOT_SET;
  MBTuplePhiId=FIELD_NOT_SET;
  MBTupleZId=FIELD_NOT_SET;
  MBTupleXId=FIELD_NOT_SET;
  MBTupleYId=FIELD_NOT_SET;
  MBTupleWorldXId=FIELD_NOT_SET;
  MBTupleWorldYId=FIELD_NOT_SET;
  MBTupleWorldZId=FIELD_NOT_SET;
  MBTupleWorldThetaId=FIELD_NOT_SET;
  MBTupleWorldPhiId=FIELD_NOT_SET;
  MBTupleFIRadiusId=FIELD_NOT_SET;
  MBTupleFIPhiId=FIELD_NOT_SET;
  MBTupleFIZId=FIELD_NOT_SET;
  MBTupleFIXId=FIELD_NOT_SET;
  MBTupleFIYId=FIELD_NOT_SET;
  MBTupleFIWorldXId=FIELD_NOT_SET;
  MBTupleFIWorldYId=FIELD_NOT_SET;
  MBTupleFIWorldZId=FIELD_NOT_SET;
  MBTupleFIWorldThetaId=FIELD_NOT_SET;
  MBTupleFIWorldPhiId=FIELD_NOT_SET;
}

//Destructor 
MiniBallHistoManager::~MiniBallHistoManager()
{
  if(HistoGenerator)
    HistoGenerator->Close();
}


//////////////////////////////////////////////////////////////////////////////
// Creating histogramms
//////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::BookHisto(G4String name, G4int cryo_id,
		G4int det_nb, G4int seg_id)
{
	//converting from relative id to absolut id
	if(cryo_id >= (int)cryostats.size())
		return;

	if(det_nb >= (int)cryostats[cryo_id]->detectorIDs.size())
		return;

	G4int det_id=cryostats[cryo_id]->detectorIDs[det_nb];

	BookHisto(name, det_id, seg_id);
}

void MiniBallHistoManager::BookHisto(G4int cryo_id, G4int det_nb, G4int seg_id)
{
	//converting from relative id to absolut id
	if(cryo_id >= (int)cryostats.size())
		return;

	if(det_nb >= (int)cryostats[cryo_id]->detectorIDs.size())
		return;

	G4int det_id=cryostats[cryo_id]->detectorIDs[det_nb];

	G4String name=
		cryostats[cryo_id]->CryostatName+"_"+detectors[det_id]->DetectorName;

	if(seg_id == CORE)
		name+="_Core";
	else if(seg_id < NB_OF_SEGMENTS)
	{
		std::ostringstream id_strm;
		id_strm.clear();
		id_strm<<seg_id;
		name+="_Seg_";
		name+=id_strm.str();
	}

	BookHisto(name, det_id, seg_id);
}

void MiniBallHistoManager::BookHisto(G4String name, G4int det_id, G4int seg_id)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	if(seg_id == CORE) 
	{
		// check if histogram allready exists
		if(detectors[det_id]->CoreHisto)
			return;

		if(detectors[det_id]->CoreEnNBin>0)
		{
			if(!HistoGenerator)
				return;

			detectors[det_id]->CoreHisto
				= HistoGenerator->NewHistogram1D(name,
						detectors[det_id]->CoreEnNBin,
						detectors[det_id]->CoreEnLow,
						detectors[det_id]->CoreEnHigh); 
			if(detectors[det_id]->CoreHisto)
				G4cout << "creating MiniBall histogram " << name << G4endl;
		}
	} 
	else if(seg_id < NB_OF_SEGMENTS) 
	{
		// check if histogram allready exists
		if(detectors[det_id]->SegmentHisto[seg_id])
			return;      

		if(detectors[det_id]->SegmentEnNBin[seg_id]>0) {
			if(!HistoGenerator)
				return;

			detectors[det_id]->SegmentHisto[seg_id]
				= HistoGenerator->NewHistogram1D(name,
						detectors[det_id]->SegmentEnNBin[seg_id],
						detectors[det_id]->SegmentEnLow[seg_id],
						detectors[det_id]->SegmentEnHigh[seg_id]);
			if(detectors[det_id]->SegmentHisto[seg_id])
				G4cout << "creating MiniBall histogram " << name << G4endl;
		}
	}
}

void MiniBallHistoManager::BookHisto(G4int det_id, G4int seg_id)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	G4String name = detectors[det_id]->DetectorName;

	if(seg_id == CORE)
		name+="_Core";
	else if(seg_id < NB_OF_SEGMENTS)
	{
		std::ostringstream id_strm;
		id_strm.clear();
		id_strm<<seg_id;
		name+="_Seg_";
		name+=id_strm.str();
	}

	BookHisto(name, det_id, seg_id);
}

void MiniBallHistoManager::BookHisto(G4String name, G4int seg_id)
{
	BookHisto(name, 0, seg_id);
}

void MiniBallHistoManager::BookHisto(G4int seg_id)
{
	//check wether at least one detector exists
	if(detectors.empty())
		return;

	G4String name="";

	if(seg_id == CORE)
		name="Core";
	else if(seg_id < NB_OF_SEGMENTS)
	{
		std::ostringstream id_strm;
		id_strm.clear();
		id_strm<<seg_id;
		name="Seg_";
		name+=id_strm.str();
	}

	BookHisto(name, 0, seg_id);
}

//////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::BookAddbackCryostat(G4String name, G4int cryo_id)
{
	// check whether cryostat exists
	if(cryo_id >= (int)cryostats.size())
		return;

	// check if histogram allready exists
	if(cryostats[cryo_id]->CryostatHisto)
		return;

	if(cryostats[cryo_id]->CryostatEnNBin>0) {
		if(!HistoGenerator)
			return;

		cryostats[cryo_id]->CryostatHisto
			= HistoGenerator->NewHistogram1D(name+"Addback",
					cryostats[cryo_id]->CryostatEnNBin,
					cryostats[cryo_id]->CryostatEnLow,
					cryostats[cryo_id]->CryostatEnHigh); 
		if(cryostats[cryo_id]->CryostatHisto)
			G4cout << "creating MiniBall histogram " << name << G4endl;
	}
}

void MiniBallHistoManager::BookAddbackCryostat(G4int cryo_id)
{
	// check whether cryostat exists
	if(cryo_id >= (int)cryostats.size())
		return;

	BookAddbackCryostat((cryostats[cryo_id]->CryostatName), cryo_id);
}

//////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::BookAllInDetector(G4String name, G4int cryo_id,
		G4int det_nb)
{
	//converting from relative id to absolut id
	if(cryo_id >= (int)cryostats.size())
		return;

	if(det_nb >= (int)cryostats[cryo_id]->detectorIDs.size())
		return;

	G4int det_id=cryostats[cryo_id]->detectorIDs[det_nb];

	BookAllInDetector(name, det_id);
}

void MiniBallHistoManager::BookAllInDetector(G4int cryo_id, G4int det_nb)
{
	//converting from relative id to absolut id
	if(cryo_id >= (int)cryostats.size())
		return;

	if(det_nb >= (int)cryostats[cryo_id]->detectorIDs.size())
		return;

	G4int det_id=cryostats[cryo_id]->detectorIDs[det_nb];

	G4String name =
		cryostats[cryo_id]->CryostatName+"_"+detectors[det_id]->DetectorName;

	BookAllInDetector(name, det_id);
}

void MiniBallHistoManager::BookAllInDetector(G4String name, G4int det_id)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	BookHisto(name + "_Core", det_id, CORE);
	for(G4int seg_id=0; seg_id<NB_OF_SEGMENTS; seg_id++)
	{
		std::ostringstream id_strm;
		id_strm.clear();
		id_strm<<seg_id;

		BookHisto(name+"_Seg_"+id_strm.str() ,det_id, seg_id);
	}
}

void MiniBallHistoManager::BookAllInDetector(G4int det_id)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	BookAllInDetector(detectors[det_id]->DetectorName, det_id);
}

void MiniBallHistoManager::BookAllInCryostat(G4String name, G4int cryo_id)
{
	// check whether cryostat exists
	if(cryo_id >= (int)cryostats.size())
		return;

	G4int nb_of_det = cryostats[cryo_id]->detectorIDs.size();
	for(G4int det_nb=0; det_nb<nb_of_det; det_nb++)
	{
		G4String det_name;
		if(nb_of_det < 26)
			det_name = G4String(det_nb+'A');
		else
		{
			std::ostringstream det_strm;
			det_strm.clear();
			det_strm << det_nb;
			det_name = det_strm.str();
		}
		BookAllInDetector(name + "_Det" + det_name, cryo_id, det_nb);
	}
}

void MiniBallHistoManager::BookAllInCryostat(G4int cryo_id)
{
	BookAllInCryostat((cryostats[cryo_id]->CryostatName), cryo_id);
}


void MiniBallHistoManager::BookAll(G4String name)
{
	for(size_t cryo_id=0; cryo_id<cryostats.size();cryo_id++)
	{
		BookAddbackCryostat(name,cryo_id);
		BookAllInCryostat(name, cryo_id);
	}
}

void MiniBallHistoManager::BookAll()
{
	for(size_t cryo_id=0; cryo_id < cryostats.size(); cryo_id++)
	{
		BookAddbackCryostat(cryo_id);
		BookAllInCryostat(cryo_id);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Creating Tuple-Fields
//////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::AddFieldSegmentEn(G4String name, G4int cryo_id,
		G4int det_nb, G4int seg_id)
{
	//converting from relative det_nb to absolut det_id
	if(cryo_id >= (int)cryostats.size())
		return;

	if(det_nb >= (int)cryostats[cryo_id]->detectorIDs.size())
		return;

	G4int det_id=cryostats[cryo_id]->detectorIDs[det_nb];

	AddFieldSegmentEn(name, det_id, seg_id);
}

void MiniBallHistoManager::AddFieldSegmentEn(G4int cryo_id, G4int det_nb,
		G4int seg_id)
{
	//converting from relative det_id to absolut det_id
	if(cryo_id >= (int)cryostats.size())
		return;

	if(det_nb >= (int)cryostats[cryo_id]->detectorIDs.size())
		return;

	G4int det_id=cryostats[cryo_id]->detectorIDs[det_nb];

	G4String name=
		cryostats[cryo_id]->CryostatName+"_"+detectors[det_id]->DetectorName;

	if(seg_id == CORE)
		name+="_Core";
	else if(seg_id < NB_OF_SEGMENTS)
	{
		std::ostringstream id_strm;
		id_strm.clear();
		id_strm<<seg_id;
		name+="_Seg_";
		name+=id_strm.str();
	}

	AddFieldSegmentEn(name, det_id, seg_id);
}

void MiniBallHistoManager::AddFieldSegmentEn(G4String name,
		G4int det_id, G4int seg_id)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	if(!MBTuple && (seg_id == CORE || seg_id<NB_OF_SEGMENTS) )
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	if(seg_id == CORE) 
	{
		if(detectors[det_id]->MBTupleCoreEnId != FIELD_NOT_SET)
			std::cout << "Warning: CORE of detector ID " << det_id
				<< " already assigned" << std::endl;
		else
		{
			detectors[det_id]->MBTupleCoreEnId = MBTuple->AddColumn(name);
			AnyMBTupleCoreEn = true;
		}
	} 
	else if(seg_id<NB_OF_SEGMENTS)
	{
		if(detectors[det_id]->MBTupleSegmentEnId[seg_id] != FIELD_NOT_SET)
			std::cout << "Warning: segment " << seg_id << " of detector ID " << det_id
				<< " already assigned" << std::endl;
		else
		{
			detectors[det_id]->MBTupleSegmentEnId[seg_id]
				= MBTuple->AddColumn(name);
			AnyMBTupleSegmentEn = true;
		}
	}
}

void MiniBallHistoManager::AddFieldSegmentEn(G4int det_id, G4int seg_id)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	G4String name=detectors[det_id]->DetectorName;

	if(seg_id == CORE)
		name+="_Core";
	else if(seg_id < NB_OF_SEGMENTS)
	{
		std::ostringstream id_strm;
		id_strm.clear();
		id_strm<<seg_id;
		name+="_Seg_";
		name+=id_strm.str();
	}

	AddFieldSegmentEn(name, det_id, seg_id);
}

void MiniBallHistoManager::AddFieldSegmentEn(G4String name, G4int seg_id)
{
	AddFieldSegmentEn(name,0,seg_id);
}

void MiniBallHistoManager::AddFieldSegmentEn(G4int seg_id)
{
	//check if at least one detector exists
	if(detectors.empty())
		return;

	AddFieldSegmentEn((detectors[0])->DetectorName,0,seg_id);
}

void MiniBallHistoManager::AddFieldAddbackCryostatEn(G4String name,
		G4int cryo_id)
{
	// check whether cryostat exists
	if(cryo_id >= (int)cryostats.size())
		return;

	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	if(cryostats[cryo_id]->MBTupleCryostatEnId != FIELD_NOT_SET)
		std::cout << "Warning: addback energy in cryostat "
			<< cryo_id << " already assigned" << std::endl;
	else
	{
		cryostats[cryo_id]->MBTupleCryostatEnId = MBTuple->AddColumn(name);
		AnyMBTupleCryostatEn = true;
	}
}

void MiniBallHistoManager::AddFieldAddbackCryostatEn(G4int cryo_id)
{
	// check whether cryostat exists
	if(cryo_id >= (int)cryostats.size())
		return;

	AddFieldAddbackCryostatEn((cryostats[cryo_id]->CryostatName), cryo_id);
}

///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::AddFieldAllInDetectorEn(G4String name, 
		G4int cryo_id, G4int det_nb)
{
	//converting from relative det_id to absolut det_id
	if(cryo_id >= (int)cryostats.size())
		return;

	if(det_nb >= (int)cryostats[cryo_id]->detectorIDs.size())
		return;

	G4int det_id=cryostats[cryo_id]->detectorIDs[det_nb];

	AddFieldAllInDetectorEn(name, det_id);
}

void MiniBallHistoManager::AddFieldAllInDetectorEn(G4int cryo_id, G4int det_nb)
{
	//converting from relative det_id to absolut det_id
	if(cryo_id >= (int)cryostats.size())
		return;

	if(det_nb >= (int)cryostats[cryo_id]->detectorIDs.size())
		return;

	G4int det_id=cryostats[cryo_id]->detectorIDs[det_nb];

	G4String name=
		cryostats[cryo_id]->CryostatName+"_"+detectors[det_id]->DetectorName;

	AddFieldAllInDetectorEn(name, det_id);
}

void MiniBallHistoManager::AddFieldAllInDetectorEn(G4String name, G4int det_id)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	detectors[det_id]->MBTupleCoreEnId
		= MBTuple->AddColumn(name + "_Core");
	for(G4int seg_id=0; seg_id<NB_OF_SEGMENTS; seg_id++)
	{
		std::ostringstream id_strm;
		id_strm.clear();
		id_strm<<seg_id;

		detectors[det_id]->MBTupleSegmentEnId[seg_id]
			= MBTuple->AddColumn(name + "_Seg_" + id_strm.str());
	}

	AnyMBTupleCoreEn = true;
	AnyMBTupleSegmentEn = true;
}

void MiniBallHistoManager::AddFieldAllInDetectorEn(G4int det_id)
{
	//check if detector exists
	if(det_id >= (int)detectors.size())
		return;

	AddFieldAllInDetectorEn((detectors[det_id])->DetectorName, det_id);
}

void MiniBallHistoManager::AddFieldAllInCryostatEn(G4String name,
		G4int cryo_id)
{
	// check whether cryostat exists
	if(cryo_id >= (int)cryostats.size())
		return;

	AddFieldAddbackCryostatEn(name + "_Addback", cryo_id);
	G4int nb_of_det = cryostats[cryo_id]->detectorIDs.size();
	for(G4int det_nb=0; det_nb<nb_of_det; det_nb++)
	{
		G4String det_name;
		if(nb_of_det < 26)
			det_name = G4String(det_nb+'A');
		else 
		{
			std::ostringstream det_strm;
			det_strm.clear();
			det_strm << det_nb;
			det_name=det_strm.str();
		}
		AddFieldAllInDetectorEn(name + "_Det" + det_name, cryo_id, det_nb);
	}
}

void MiniBallHistoManager::AddFieldAllInCryostatEn(G4int cryo_id)
{
	// check whether cryostat exists
	if(cryo_id >= (int)cryostats.size())
		return;

	AddFieldAllInCryostatEn((cryostats[cryo_id]->CryostatName), cryo_id);
}

void MiniBallHistoManager::AddFieldAllEn(G4String name)
{
	for(size_t cryo_id=0;cryo_id<cryostats.size();cryo_id++)
	{
		AddFieldAllInCryostatEn(name, cryo_id);
	}
}

void MiniBallHistoManager::AddFieldAllEn()
{
	for(size_t cryo_id=0;cryo_id<cryostats.size();cryo_id++)
	{
		AddFieldAllInCryostatEn((cryostats[cryo_id])->CryostatName, cryo_id);
	}
}

///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::AddFieldSegmentNb(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleSegmentNbId = MBTuple->AddColumn(name, -1);
}

void MiniBallHistoManager::AddFieldSegmentCode(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleSegmentCodeId = MBTuple->AddColumn(name, -1);
}

void MiniBallHistoManager::AddFieldDetectorNb(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleDetectorNbId = MBTuple->AddColumn(name, -1);
}

void MiniBallHistoManager::AddFieldCryostatNb(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleCryostatNbId = MBTuple->AddColumn(name, -1);
}

void MiniBallHistoManager::AddFieldEventNb(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleEventNbId = MBTuple->AddColumn(name, -1);
}

void MiniBallHistoManager::AddFieldTraceNb(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleTraceNbId = MBTuple->AddColumn(name, -1);
}



///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::AddFieldPositionInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	AddFieldRadiusInDet(name+"R");
	AddFieldPhiInDet(name+"Phi");
	AddFieldZInDet(name+"Z");
}

void MiniBallHistoManager::AddFieldRadiusInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleRadiusId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldPhiInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTuplePhiId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldZInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleZId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldXInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleXId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldYInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleYId = MBTuple->AddColumn(name);
}

///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::AddFieldPositionInWorld(G4String name)
{
	AddFieldWorldX(name+"X");
	AddFieldWorldY(name+"Y");
	AddFieldWorldZ(name+"Z");
}

void MiniBallHistoManager::AddFieldWorldX(G4String name)
{
	if(!MBTuple) 
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleWorldXId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldWorldY(G4String name)
{
	if(!MBTuple) 
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleWorldYId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldWorldZ(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleWorldZId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldWorldTheta(G4String name)
{
	if(!MBTuple) 
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleWorldThetaId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldWorldPhi(G4String name)
{
	if(!MBTuple) 
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleWorldPhiId = MBTuple->AddColumn(name);
}

///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::AddFieldFIPositionInDet(G4String name)
{
	AddFieldFIRadiusInDet(name+"r");
	AddFieldFIPhiInDet(name+"phi");
	AddFieldFIZInDet(name+"z");
}

void MiniBallHistoManager::AddFieldFIRadiusInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIRadiusId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldFIPhiInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIPhiId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldFIZInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIZId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldFIXInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIXId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldFIYInDet(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIYId = MBTuple->AddColumn(name);
}

///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::AddFieldFIPositionInWorld(G4String name)
{
	AddFieldFIWorldX(name+"X");
	AddFieldFIWorldY(name+"Y");
	AddFieldFIWorldZ(name+"Z");
}

void MiniBallHistoManager::AddFieldFIWorldX(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIWorldXId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldFIWorldY(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIWorldYId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldFIWorldZ(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIWorldZId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldFIWorldTheta(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIWorldThetaId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddFieldFIWorldPhi(G4String name)
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	MBTupleFIWorldPhiId = MBTuple->AddColumn(name);
}

void MiniBallHistoManager::AddUserDefinedArray(G4String name,
		G4int nb_of_entries, G4int* segment_nb)
{
	if(!MBTuple) 
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	// create new array
	MiniBallArray* new_array = new MiniBallArray();

	for(G4int entr_nb=0; entr_nb<nb_of_entries; entr_nb++)
		new_array->segments.push_back(segment_nb[entr_nb]);

	// create array in the tuple
	G4int nb_of_fields = detectors.size()*nb_of_entries;
	new_array->MBTupleArrayID = MBTuple->AddArray(name, nb_of_fields);

	new_array->MBTupleHitSegmentID=FIELD_NOT_SET;
	new_array->MBTupleRadiusInDetID=FIELD_NOT_SET;
	new_array->MBTuplePhiInDetID=FIELD_NOT_SET;
	new_array->MBTupleZInDetID=FIELD_NOT_SET;
	new_array->MBTupleXInDetID=FIELD_NOT_SET;
	new_array->MBTupleYInDetID=FIELD_NOT_SET;

	// register new array
	arrays.push_back(new_array);
}

void MiniBallHistoManager::AddCoresArray(G4String name)
{
	G4int seg_list[] = { CORE };

	AddUserDefinedArray(name, 1, seg_list);
}

void MiniBallHistoManager::AddArrayHitSegment(G4String name)
{
	G4int last_array=arrays.size();
	if(last_array>0)
	{
		G4int nb_of_fields
			= detectors.size()*arrays[last_array-1]->segments.size();
		arrays[last_array-1]->MBTupleHitSegmentID
			= MBTuple->AddArray(name, nb_of_fields);
	}
}

void MiniBallHistoManager::AddArrayRadiusInDet(G4String name)
{
	G4int last_array=arrays.size();
	if(last_array>0)
	{
		G4int nb_of_fields
			= detectors.size()*arrays[last_array-1]->segments.size();
		arrays[last_array-1]->MBTupleRadiusInDetID
			= MBTuple->AddArray(name, nb_of_fields);
	}
}

void MiniBallHistoManager::AddArrayPhiInDet(G4String name)
{
	G4int last_array=arrays.size();
	if(last_array>0)
	{
		G4int nb_of_fields
			= detectors.size()*arrays[last_array-1]->segments.size();
		arrays[last_array-1]->MBTuplePhiInDetID
			= MBTuple->AddArray(name, nb_of_fields);
	}
}

void MiniBallHistoManager::AddArrayZInDet(G4String name)
{
	G4int last_array=arrays.size();
	if(last_array>0)
	{
		G4int nb_of_fields
			= detectors.size()*arrays[last_array-1]->segments.size();
		arrays[last_array-1]->MBTupleZInDetID
			= MBTuple->AddArray(name, nb_of_fields);
	}
}

void MiniBallHistoManager::AddArrayPositionInDet(G4String name)
{
	AddArrayRadiusInDet(name+"r");
	AddArrayPhiInDet(name+"phi");
	AddArrayZInDet(name+"z");
}

void MiniBallHistoManager::AddArrayXInDet(G4String name)
{
	G4int last_array=arrays.size();
	if(last_array>0)
	{
		G4int nb_of_fields
			= detectors.size()*arrays[last_array-1]->segments.size();
		arrays[last_array-1]->MBTupleXInDetID
			= MBTuple->AddArray(name, nb_of_fields);
	}
}

void MiniBallHistoManager::AddArrayYInDet(G4String name)
{
	G4int last_array=arrays.size();
	if(last_array>0)
	{
		G4int nb_of_fields
			= detectors.size()*arrays[last_array-1]->segments.size();
		arrays[last_array-1]->MBTupleYInDetID
			= MBTuple->AddArray(name, nb_of_fields);
	}
}

void MiniBallHistoManager::SetEnergyUnit(G4double p_energy_unit)
{
	EnergyUnit = p_energy_unit;
}

void MiniBallHistoManager::SetSigma(G4int det_id, G4int seg_id,
		G4double electronic_sigma,
		G4double statistic_epsilon)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	WriteGaussian=true;
	if(seg_id == CORE)
	{
		detectors[det_id]->CoreElectronicSigma = electronic_sigma;
		detectors[det_id]->CoreStatisticEpsilon = statistic_epsilon;
	}
	else if(seg_id < NB_OF_SEGMENTS && seg_id >= 0)
	{
		detectors[det_id]->SegmentElectronicSigma[seg_id] = electronic_sigma;
		detectors[det_id]->SegmentStatisticEpsilon[seg_id] = statistic_epsilon;
	}
	else if(seg_id == ALL_SEG)
		for(int i=0; i<NB_OF_SEGMENTS; i++)
		{
			detectors[det_id]->SegmentElectronicSigma[i] = electronic_sigma;
			detectors[det_id]->SegmentStatisticEpsilon[i] = statistic_epsilon;
		}
}  

void MiniBallHistoManager::SetAllSigma(G4int seg_id, G4double electronic_sigma,
		G4double statistic_epsilon)
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++){
		SetSigma(det_id, seg_id, electronic_sigma, statistic_epsilon);
	}
}

void MiniBallHistoManager::SetAllSigma(G4double eletronic_sigma,
		G4double statistic_epsilon)
{
	SetAllSigma(CORE, eletronic_sigma, statistic_epsilon);
	SetAllSigma(ALL_SEG, eletronic_sigma, statistic_epsilon);
}

void MiniBallHistoManager::SetThreshold(G4int det_id, G4int seg_id,
		G4double cutoff_en,
		G4double cutoff_width,
		G4double suppression_en,
		G4double suppression_slope)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	ElectronicThreshold=true;
	if(seg_id == CORE)
	{
		detectors[det_id]->CoreCutoffEn = cutoff_en;
		detectors[det_id]->CoreCutoffWidth = cutoff_width;
		detectors[det_id]->CoreSuppressionEn = suppression_en;
		detectors[det_id]->CoreSuppressionSlope = suppression_slope;
	}
	else if(seg_id < NB_OF_SEGMENTS && seg_id >= 0)
	{
		detectors[det_id]->SegmentCutoffEn[seg_id] = cutoff_en;
		detectors[det_id]->SegmentCutoffWidth[seg_id] = cutoff_width;
		detectors[det_id]->SegmentSuppressionEn[seg_id] = suppression_en;
		detectors[det_id]->SegmentSuppressionSlope[seg_id] = suppression_slope;
	}
	else if(seg_id == ALL_SEG)
		for(int i=0; i<NB_OF_SEGMENTS; i++)
		{
			detectors[det_id]->SegmentCutoffEn[i] = cutoff_en;
			detectors[det_id]->SegmentCutoffWidth[i] = cutoff_width;
			detectors[det_id]->SegmentSuppressionEn[i] = suppression_en;
			detectors[det_id]->SegmentSuppressionSlope[i] = suppression_slope;
		}
}  

void MiniBallHistoManager::SetAllThreshold(G4int seg_id, G4double cutoff_en,
		G4double cutoff_width,
		G4double suppression_en,
		G4double suppression_slope)
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++)
		SetThreshold(det_id, seg_id, cutoff_en, cutoff_width,
				suppression_en, suppression_slope);
}

void MiniBallHistoManager::SetAllThreshold(G4double cutoff_en,
		G4double cutoff_width,
		G4double suppression_en,
		G4double suppression_slope)
{
	SetAllThreshold(CORE, cutoff_en, cutoff_width,
			suppression_en, suppression_slope);
	SetAllThreshold(ALL_SEG, cutoff_en, cutoff_width,
			suppression_en, suppression_slope);
}

///////////////////////////////////////////////////////////////////////////////
// intern organization
//////////////////////////////////////////////////////////////////////////////

G4int MiniBallHistoManager::AddDetector(G4String name)
{
	MiniBallDetectorHisto* new_detector
		= new MiniBallDetectorHisto(name);

	// initialize detector values
	new_detector->CoreElectronicSigma=0.;
	new_detector->CoreStatisticEpsilon=0.;
	new_detector->CoreCutoffEn=0.;
	new_detector->CoreCutoffWidth=-1.;
	new_detector->CoreSuppressionEn=0.;
	new_detector->CoreSuppressionSlope=-1.;
	for(int i=0; i<NB_OF_SEGMENTS; i++)
	{
		new_detector->SegmentElectronicSigma[i]=0.;
		new_detector->SegmentStatisticEpsilon[i]=0.;
		new_detector->SegmentCutoffEn[i]=0.;
		new_detector->SegmentCutoffWidth[i]=-1.;
		new_detector->SegmentSuppressionEn[i]=0.;
		new_detector->SegmentSuppressionSlope[i]=-1.;
	}

	G4int det_id = detectors.size();
	detectors.push_back(new_detector);

	return(det_id);
}

G4int MiniBallHistoManager::AddDetector()
{
	G4int det_id = detectors.size();
	std::ostringstream id_strm;
	id_strm.clear();
	id_strm << det_id;
	return( AddDetector("det" + id_strm.str()) );
}

G4int MiniBallHistoManager::AddCryostat(G4String name)
{
	MiniBallCryostatHisto* new_cryostat
		=new MiniBallCryostatHisto(name);


	size_t cryo_id = cryostats.size();
	cryostats.push_back(new_cryostat);

	return(cryo_id);
}

G4int MiniBallHistoManager::AddCryostat()
{
	size_t cryo_id = cryostats.size();
	std::ostringstream id_strm;
	id_strm.clear();
	id_strm << cryo_id;
	return( AddCryostat("cryo" + id_strm.str()) );
}

void MiniBallHistoManager::AddToCryostat(G4int cryo_id, G4int det_id)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	// check whether cryostat exists
	if(cryo_id >= (int)cryostats.size())
		return;  

	//if cryo_id and det_id are valid, 
	//then the new detector "det_id" is added to the cryostat "cryo_id"
	cryostats[cryo_id]->detectorIDs.push_back(det_id);
	detectors[det_id]->cryostatID=cryo_id;
}

///////////////////////////////////////////////////////////////////////////////

G4int MiniBallHistoManager::GetNbOfDetectors()
{
	return(detectors.size());
}

G4int MiniBallHistoManager::GetNbOfCryostats()
{
	return(cryostats.size());
}

///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::SetMBHitsCollectionID(G4int det_id, G4int HCID)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size() || det_id < 0)
		return;
	else
		detectors[det_id]->DetectorHitsCollectionID=HCID;
}

G4int MiniBallHistoManager::GetMBHitsCollectionID(G4int det_id)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size() || det_id < 0)
		return(-1);
	else
		return(detectors[det_id]->DetectorHitsCollectionID);
}

///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::SetMBEnNBin(G4int det_id,
		G4int seg_id, G4int NBin)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	if(seg_id == CORE) detectors[det_id]->CoreEnNBin = NBin;
	else if(seg_id < NB_OF_SEGMENTS)
		detectors[det_id]->SegmentEnNBin[seg_id] = NBin;
}

void MiniBallHistoManager::SetMBEnLow(G4int det_id,
		G4int seg_id, G4double Low)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;
	if(seg_id == CORE) detectors[det_id]->CoreEnLow = Low;
	else if(seg_id < NB_OF_SEGMENTS && seg_id >= 0)
		detectors[det_id]->SegmentEnLow[seg_id] = Low;
}

void MiniBallHistoManager::SetMBEnHigh(G4int det_id,
		G4int seg_id, G4double High)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	if(seg_id == CORE) detectors[det_id]->CoreEnHigh = High;
	else if(seg_id < NB_OF_SEGMENTS && seg_id >= 0)
		detectors[det_id]->SegmentEnHigh[seg_id] = High;
}

void MiniBallHistoManager::SetMBEnNBin(G4int det_id, G4int NBin)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	SetMBEnNBin(det_id, CORE, NBin);
	for(G4int seg_id=0; seg_id<NB_OF_SEGMENTS; seg_id++)
		SetMBEnNBin(det_id, seg_id, NBin);
}

void MiniBallHistoManager::SetMBEnLow(G4int det_id, G4double Low)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	SetMBEnLow(det_id, CORE, Low);
	for(G4int seg_id=0; seg_id<NB_OF_SEGMENTS; seg_id++)
		SetMBEnLow(det_id, seg_id, Low);
}

void MiniBallHistoManager::SetMBEnHigh(G4int det_id, G4double High)
{
	// check whether detector exists
	if(det_id >= (int)detectors.size())
		return;

	SetMBEnHigh(det_id, CORE, High);
	for(G4int seg_id=0; seg_id<NB_OF_SEGMENTS; seg_id++)
		SetMBEnHigh(det_id, seg_id, High);
}

void MiniBallHistoManager::SetMBCryostatEnNBin(G4int cryo_id, G4int NBin)
{
	// check whether detector exists
	if(cryo_id >= (int)cryostats.size())
		return;

	cryostats[cryo_id]->CryostatEnNBin = NBin;
}

void MiniBallHistoManager::SetMBCryostatEnLow(G4int cryo_id, G4double Low)
{
	// check whether detector exists
	if(cryo_id >= (int)cryostats.size())
		return;

	cryostats[cryo_id]->CryostatEnLow = Low;
}

void MiniBallHistoManager::SetMBCryostatEnHigh(G4int cryo_id, G4double High)
{
	// check whether cryostat exists
	if(cryo_id >= (int)cryostats.size())
		return;

	cryostats[cryo_id]->CryostatEnHigh = High;
}

void MiniBallHistoManager::SetMBEnNBin(G4int NBin)
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++)
		SetMBEnNBin(det_id, NBin);
	for(size_t cryo_id=0; cryo_id<cryostats.size(); cryo_id++)
		SetMBCryostatEnNBin(cryo_id, NBin);
}

void MiniBallHistoManager::SetMBEnLow(G4double Low)
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++)
		SetMBEnLow(det_id, Low);
	for(size_t cryo_id=0; cryo_id<cryostats.size(); cryo_id++)
		SetMBCryostatEnLow(cryo_id, Low);
}   

void MiniBallHistoManager::SetMBEnHigh(G4double High)
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++)
		SetMBEnHigh(det_id, High);
	for(size_t cryo_id=0; cryo_id<cryostats.size(); cryo_id++)
		SetMBCryostatEnHigh(cryo_id, High);
}   

///////////////////////////////////////////////////////////////////////////////

MiniBallTuple* MiniBallHistoManager::GetMiniBallTuple()
{
	if(!MBTuple)
		MBTuple = HistoGenerator->NewTuple(MBTupleName);

	return(MBTuple);
}

void MiniBallHistoManager::SetMiniBallTuple(MiniBallTuple* tup)
{
	MBTuple=tup;
	MBTupleManagement=EXTERN_MANAGEMENT;
}

///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::SetMBInteraction(size_t det_id, size_t seg_id,
		G4double en,
		G4int event_id,
		G4int trace_id,
		G4ThreeVector positionInDet,
		G4ThreeVector positionInWorld,
		G4double time)
{
	//First Interaction in this event:
	if(!InteractionInEvent)
	{
		FirstInteractionPositionInDet.set(positionInDet.x(),
				positionInDet.y(), 
				positionInDet.z());
		FirstInteractionPositionInWorld.set(positionInWorld.x(),
				positionInWorld.y(),
				positionInWorld.z());
	}

	Interaction=true;
	InteractionInEvent=true;

	//This single interaction:
	InteractionEn=en;
	InteractionDetectorNb=det_id;
	InteractionEventNb=event_id;
	InteractionTraceNb=trace_id;
	InteractionSegmentNb=seg_id;
	InteractionPositionInDet.set(positionInDet.x(),
			positionInDet.y(),
			positionInDet.z());
	InteractionPositionInWorld.set(positionInWorld.x(), 
			positionInWorld.y(),
			positionInWorld.z());


	//main interaction in this event
	if(InteractionEn>MainInteractionEn) 
	{
		MainInteractionEn=en;
		MainInteractionDetectorNb=det_id;
		MainInteractionSegmentNb=seg_id;
		MainInteractionPositionInDet.set(positionInDet.x(),
				positionInDet.y(),
				positionInDet.z());
		MainInteractionPositionInWorld.set(positionInWorld.x(),
				positionInWorld.y(),
				positionInWorld.z());
		if(det_id < detectors.size())
		{ 
			InteractionCryostatNb=detectors[det_id]->cryostatID;
			MainInteractionCryostatNb=InteractionCryostatNb;
		}   
	}

	//main interaction for the hit detector
	if(InteractionEn>detectors[det_id]->MIEnergy)
	{
		detectors[det_id]->MIEnergy=en;
		detectors[det_id]->MIPositionInDet.set(positionInDet.x(),
				positionInDet.y(),
				positionInDet.z());
		detectors[det_id]->TRInteractionPosition.set(positionInWorld.x(),
				positionInWorld.y(),
				positionInWorld.z());


		detectors[det_id]->MIHitSegment=seg_id;
		detectors[det_id]->time = time;
	}

	//energy deposit in this interaction
	if(seg_id<(size_t)NB_OF_SEGMENTS)
	{
		// fill core
		detectors[det_id]->CoreEn=en;
		detectors[det_id]->CoreEnAccu+=en;
		// fill segment
		detectors[det_id]->SegmentEn[seg_id]=en;
		detectors[det_id]->SegmentEnAccu[seg_id]+=en;
	}
}

///////////////////////////////////////////////////////////////////////////////

void MiniBallHistoManager::StartOfEvent()
{
	InitMBAccumulation();
	InteractionInEvent=false;

	if(MBTuple)
		MBTuple->Clear();
	if(MBTupleType==EVENT) 
		WroteMBTuple = false;
}

void MiniBallHistoManager::EndOfEvent()
{
	if(WriteGaussian)
		AddGaussian();
	if(ElectronicThreshold)
		SetCutoff();
	FillMBHisto();
	if(MBTupleType==EVENT && (InteractionInEvent || WroteMBTuple))
		FillMBTuple();
}

void MiniBallHistoManager::StartOfInteraction()
{
	InitMBInteraction();
	Interaction=false;

	if(MBTupleType==INTERACTION) 
		WroteMBTuple = false;
}

void MiniBallHistoManager::EndOfInteraction()
{
	if(MBTupleType==INTERACTION && (Interaction || WroteMBTuple))
		FillMBInteractionTuple();
}

///////////////////////////////////////////////////////////////////////////////

//private:

// resets all entries for accumulation of energy during an event
void MiniBallHistoManager::InitMBAccumulation()
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++)
	{
		detectors[det_id]->CoreEnAccu=0;
		detectors[det_id]->MIEnergy=0;
		detectors[det_id]->MIHitSegment=-1;
		for(G4int segment_id=0; segment_id<NB_OF_SEGMENTS; segment_id++)
			detectors[det_id]->SegmentEnAccu[segment_id]=0;
	}

	MainInteractionEn=0;
	MainInteractionDetectorNb=-1;
	MainInteractionSegmentNb=-1;
	MainInteractionCryostatNb=-1;

	MainInteractionPositionInDet.set(0., 0., 0.);
	MainInteractionPositionInWorld.set(0., 0., 0.);
	FirstInteractionPositionInDet.set(0., 0., 0.);
	FirstInteractionPositionInWorld.set(0., 0., 0.);
}

//resets all entries 
void MiniBallHistoManager::InitMBInteraction()
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++)
	{
		detectors[det_id]->CoreEn = 0;
		for(G4int segment_id=0; segment_id<NB_OF_SEGMENTS; segment_id++)
			detectors[det_id]->SegmentEn[segment_id] = 0;
	}

	InteractionEn=0;
	InteractionDetectorNb=-1;
	InteractionSegmentNb=-1;
	InteractionCryostatNb=-1;
	InteractionEventNb=-1;
	InteractionTraceNb=-1;

	InteractionPositionInDet.set(0., 0., 0.);
	InteractionPositionInWorld.set(0., 0., 0.);
}

///////////////////////////////////////////////////////////////////////////////

// add gaussian
void MiniBallHistoManager::AddGaussian()
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++)
	{
		G4double el_s, st_e;
		el_s = detectors[det_id]->CoreElectronicSigma;
		st_e = detectors[det_id]->CoreStatisticEpsilon;
		if((el_s>0 || st_e>0) && detectors[det_id]->CoreEnAccu>0)
			detectors[det_id]->CoreEnAccu
				= RandSigma(el_s, st_e, detectors[det_id]->CoreEnAccu);
		for(int seg_id=0; seg_id<NB_OF_SEGMENTS; seg_id++)
		{
			el_s = detectors[det_id]->SegmentElectronicSigma[seg_id];
			st_e = detectors[det_id]->SegmentStatisticEpsilon[seg_id];
			if((el_s>0 || st_e>0) && detectors[det_id]->SegmentEnAccu[seg_id]>0)
				detectors[det_id]->SegmentEnAccu[seg_id]
					= RandSigma(el_s,st_e,detectors[det_id]->SegmentEnAccu[seg_id]);
		}
	}
}	

// see Glenn F. Knoll "Radiation Detection and Measurement", ch. 12, III, B 
G4double MiniBallHistoManager::RandSigma(G4double electronic_sigma,
		G4double statistic_epsilon,
		G4double en)
{
	G4double W_E_2 = electronic_sigma*electronic_sigma;
	G4double W_D_2 = statistic_epsilon*en;

	G4double W_T = sqrt(W_E_2+W_D_2);

	return(en+CLHEP::RandGauss::shoot(0., W_T));
}

// cutoff
void MiniBallHistoManager::SetCutoff()
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++)
	{
		G4double cutoff_en, cutoff_width, sup_en, sup_slope, en;
		en = detectors[det_id]->CoreEnAccu;
		cutoff_en = detectors[det_id]->CoreCutoffEn;
		cutoff_width = detectors[det_id]->CoreCutoffWidth;
		sup_en = detectors[det_id]->CoreSuppressionEn;
		sup_slope = detectors[det_id]->CoreSuppressionSlope;
		if(cutoff_width>0 && en>0)
		{
			if(RandCutoff(en, cutoff_en, cutoff_width, sup_en, sup_slope))
				detectors[det_id]->CoreEnAccu = 0;
		}
		for(int seg_id=0; seg_id<NB_OF_SEGMENTS; seg_id++)
		{
			en = detectors[det_id]->SegmentEnAccu[seg_id];
			cutoff_en = detectors[det_id]->SegmentCutoffEn[seg_id];
			cutoff_width = detectors[det_id]->SegmentCutoffWidth[seg_id];
			sup_en = detectors[det_id]->SegmentSuppressionEn[seg_id];
			sup_slope = detectors[det_id]->SegmentSuppressionSlope[seg_id];
			if(cutoff_width>0 && en>0)
			{
				if(RandCutoff(en, cutoff_en, cutoff_width, sup_en, sup_slope))
					detectors[det_id]->SegmentEnAccu[seg_id] = 0;
			}
		}
	}
}	

G4bool MiniBallHistoManager::RandCutoff(G4double en, G4double cutoff_en,
		G4double cutoff_width,
		G4double sup_en, G4double sup_slope)
{
	if(sup_slope<=0 && en > cutoff_en+10*cutoff_width)
		return(false);
	if(sup_slope>0 && en > sup_en+10*cutoff_width)
		return(false);

	G4double rand = CLHEP::RandFlat::shoot(0., 1.);

	G4double value;
	if(sup_slope<=0)
		value = 0.5*(erf((en-cutoff_en)/cutoff_width)+1);
	else
	{
		G4double lin = sup_slope*(en-sup_en+1./3.*cutoff_width)+1;
		G4double erf1
			= 0.25*(erf((en-cutoff_en)/cutoff_width)+1)
			*(erf((sup_en-en)/cutoff_width)+1);
		G4double erf2 = 0.5*(erf((en-sup_en)/cutoff_width)+1);
		value = lin*erf1+erf2;
	}

	if(rand < value)
		return(false);
	else
		return(true);  
}

// fill the accumulated energies in histogram
void MiniBallHistoManager::FillMBHisto()
{
	for(size_t det_id=0; det_id<detectors.size(); det_id++)
	{
		if(detectors[det_id]->CoreEnAccu > 0)
			if(detectors[det_id]->CoreHisto)
				detectors[det_id]->CoreHisto->Fill(
						detectors[det_id]->CoreEnAccu/EnergyUnit);
		for(size_t segment_id=0; segment_id<(size_t)NB_OF_SEGMENTS; segment_id++)
			if(detectors[det_id]->SegmentEnAccu[segment_id] > 0)
				if(detectors[det_id]->SegmentHisto[segment_id])
					detectors[det_id]->SegmentHisto[segment_id]->Fill(
							detectors[det_id]->SegmentEnAccu[segment_id]/EnergyUnit);
	}

	for(size_t cryo_id=0; cryo_id<cryostats.size(); cryo_id++)
		if(cryostats[cryo_id]->CryostatHisto)
		{
			size_t nb_of_dets=cryostats[cryo_id]->detectorIDs.size();
			G4double cryo_en=0;
			for(size_t det_id=0; det_id<nb_of_dets; det_id++)
				if(detectors[det_id]->CoreEnAccu > 0)
					cryo_en += detectors[det_id]->CoreEnAccu;
			if(cryo_en>0)
				cryostats[cryo_id]->CryostatHisto->Fill(cryo_en/EnergyUnit);
		}
}  

void MiniBallHistoManager::FillMBTuple()
{
	if(!MBTuple)
		return;

	// fill detector arrays
	if(arrays.size() > 0)
	{
		for(size_t det_id=0; det_id<detectors.size(); det_id++)
		{
			// fill arrays, for each detector
			for(size_t array_nb=0; array_nb<arrays.size(); array_nb++)
			{
				size_t nb_of_entries = arrays[array_nb]->segments.size();
				G4int array_id = arrays[array_nb]->MBTupleArrayID;

				for(size_t entry_nb=0; entry_nb<nb_of_entries; entry_nb++)
				{
					size_t seg_id = arrays[array_nb]->segments[entry_nb];

					G4double value=0;
					if(seg_id == CORE)
						value = detectors[det_id]->CoreEnAccu;
					else if(seg_id<(size_t)NB_OF_SEGMENTS)
						value = detectors[det_id]->SegmentEnAccu[seg_id];

					G4int index = entry_nb+det_id*nb_of_entries;
					MBTuple->FillArray(array_id, index, value/EnergyUnit);
					if(seg_id == CORE)
					{
						if(arrays[array_nb]->MBTupleHitSegmentID>0)
							MBTuple->FillArray(
									arrays[array_nb]->MBTupleHitSegmentID, 
									index, detectors[det_id]->MIHitSegment);

						G4ThreeVector hit_position
							= detectors[det_id]->MIPositionInDet;

						if(arrays[array_nb]->MBTupleRadiusInDetID>0)
							MBTuple->FillArray(
									arrays[array_nb]->MBTupleRadiusInDetID, 
									index, hit_position.rho());

						if(arrays[array_nb]->MBTuplePhiInDetID>0)
							MBTuple->FillArray(
									arrays[array_nb]->MBTuplePhiInDetID, 
									index, hit_position.phi());

						if(arrays[array_nb]->MBTupleZInDetID>0)
							MBTuple->FillArray(
									arrays[array_nb]->MBTupleZInDetID, 
									index, hit_position.z());

						if(arrays[array_nb]->MBTupleXInDetID>0)
							MBTuple->FillArray(
									arrays[array_nb]->MBTupleXInDetID, 
									index, hit_position.x());

						if(arrays[array_nb]->MBTupleYInDetID>0)
							MBTuple->FillArray(
									arrays[array_nb]->MBTupleYInDetID, 
									index, hit_position.y());		      
					}
				}
			}
		}
	}

	// fill ntuple with DetectorEn
	if(AnyMBTupleCoreEn)
	{
		for(size_t det_id=0; det_id<detectors.size(); det_id++)
		{
			G4int column_id = detectors[det_id]->MBTupleCoreEnId;
			if(column_id!=FIELD_NOT_SET)
			{
				G4double EnAccu = detectors[det_id]->CoreEnAccu;

				if(EnAccu>0)
				{
					MBTuple->FillColumn(column_id, EnAccu/EnergyUnit);
					WroteMBTuple=true;
				}
			}
		}
	}  	     
	// fill ntuple for SegmentEn
	if(AnyMBTupleSegmentEn)
	{
		for(size_t det_id=0; det_id<detectors.size(); det_id++)
		{
			for(size_t seg=0; seg<(size_t)NB_OF_SEGMENTS; seg++)
			{
				G4int column_id = detectors[det_id]->MBTupleSegmentEnId[seg];
				if(column_id!=FIELD_NOT_SET)
				{
					G4double EnAccu = detectors[det_id]->SegmentEnAccu[seg];
					if(EnAccu>0)
					{
						MBTuple->FillColumn(column_id, EnAccu/EnergyUnit);
						WroteMBTuple=true;
					}
				}
			}
		}
	}
	// fill ntuple for CryostatEn
	if(AnyMBTupleCryostatEn)
	{
		for(size_t cryo_id=0; cryo_id<cryostats.size(); cryo_id++)
		{
			G4int column_id=cryostats[cryo_id]->MBTupleCryostatEnId;
			if(column_id!=FIELD_NOT_SET)
			{
				size_t nb_of_dets=cryostats[cryo_id]->detectorIDs.size();
				G4double cryo_en=0;
				for(size_t det_nb=0; det_nb<nb_of_dets; det_nb++)
				{
					size_t det_id = cryostats[cryo_id]->detectorIDs[det_nb];
					if(detectors[det_id]->CoreEnAccu > 0)
						cryo_en += detectors[det_id]->CoreEnAccu;
				}
				if(cryo_en>0)
				{
					MBTuple->FillColumn(column_id, cryo_en/EnergyUnit);
					WroteMBTuple=true;
				}
			}
		}
	}

	// fill ntuple for main interaction
	if(MBTupleSegmentNbId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleSegmentNbId,
				(G4double)MainInteractionSegmentNb);
		WroteMBTuple=true;
	}
	if(MBTupleSegmentCodeId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		G4int segmentCode
			= MiniBallSingleDetector::GetSegmentCode(MainInteractionSegmentNb); 
		MBTuple->FillColumn(MBTupleSegmentCodeId, (G4double)segmentCode);
		WroteMBTuple=true;
	}
	if(MBTupleDetectorNbId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleDetectorNbId,
				(G4double)MainInteractionDetectorNb);
		WroteMBTuple=true;
	}
	if(MBTupleCryostatNbId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleCryostatNbId,
				(G4double)MainInteractionCryostatNb);
		WroteMBTuple=true;
	}
	if(MBTupleRadiusId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleRadiusId, MainInteractionPositionInDet.rho());
		WroteMBTuple=true;
	}
	if(MBTuplePhiId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTuplePhiId, MainInteractionPositionInDet.phi());
		WroteMBTuple=true;
	}
	if(MBTupleZId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleZId, MainInteractionPositionInDet.z());
		WroteMBTuple=true;
	}
	if(MBTupleXId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleXId, MainInteractionPositionInDet.x());
		WroteMBTuple=true;
	}
	if(MBTupleYId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleYId, MainInteractionPositionInDet.y());
		WroteMBTuple=true;
	}
	if(MBTupleWorldXId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleWorldXId, MainInteractionPositionInWorld.x());
		WroteMBTuple=true;
	}
	if(MBTupleWorldYId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleWorldYId, MainInteractionPositionInWorld.y());
		WroteMBTuple=true;
	}
	if(MBTupleWorldZId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleWorldZId, MainInteractionPositionInWorld.z());
		WroteMBTuple=true;
	}
	if(MBTupleWorldThetaId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleWorldThetaId,
				MainInteractionPositionInWorld.theta());
		WroteMBTuple=true;
	}
	if(MBTupleWorldPhiId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleWorldPhiId,
				MainInteractionPositionInWorld.phi());
		WroteMBTuple=true;
	}
	//Fill for first Interaction
	if(MBTupleFIRadiusId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleFIRadiusId,
				FirstInteractionPositionInDet.rho());
		WroteMBTuple=true;
	}
	if(MBTupleFIPhiId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleFIPhiId, FirstInteractionPositionInDet.phi());
		WroteMBTuple=true;
	}
	if(MBTupleFIZId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleFIZId, FirstInteractionPositionInDet.z());
		WroteMBTuple=true;
	}
	if(MBTupleFIXId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleFIXId, FirstInteractionPositionInDet.x());
		WroteMBTuple=true;
	}
	if(MBTupleFIYId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleFIYId, FirstInteractionPositionInDet.y());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldXId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleFIWorldXId,
				FirstInteractionPositionInWorld.x());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldYId!=FIELD_NOT_SET && InteractionInEvent) 
	{
		MBTuple->FillColumn(MBTupleFIWorldYId,
				FirstInteractionPositionInWorld.y());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldZId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIWorldZId,
				FirstInteractionPositionInWorld.z());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldThetaId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIWorldThetaId,
				FirstInteractionPositionInWorld.theta());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldPhiId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIWorldPhiId,
				FirstInteractionPositionInWorld.phi());
		WroteMBTuple=true;
	}
	//Fill EventNb
	if(MBTupleEventNbId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleEventNbId,(G4double)InteractionEventNb);
		WroteMBTuple=true;
	}

	if(MBTupleManagement==INTERN_MANAGEMENT) {
		if(WroteMBTuple){
			MBTuple->Fill();
		}
		else{
			MBTuple->Clear();
		}
	}
}

void MiniBallHistoManager::FillMBInteractionTuple()
{
	if(!MBTuple)
		return;

	// fill ntuple with DetectorEn
	if(AnyMBTupleCoreEn)
		for(size_t det_id=0; det_id<detectors.size(); det_id++)
		{
			G4int column_id = detectors[det_id]->MBTupleCoreEnId;
			if(column_id!=FIELD_NOT_SET) 
			{
				G4double En = detectors[det_id]->CoreEn;
				if(En>0)
				{
					MBTuple->FillColumn(column_id, En/EnergyUnit);
					WroteMBTuple=true;
				}
			}
		}

	// fill ntuple for SegmentEn
	if(AnyMBTupleSegmentEn)
		for(size_t det_id=0; det_id<detectors.size(); det_id++)
			for(size_t seg=0; seg<(size_t)NB_OF_SEGMENTS; seg++) 
			{
				G4int column_id = detectors[det_id]->MBTupleSegmentEnId[seg];
				if(column_id!=FIELD_NOT_SET) 
				{
					G4double En = detectors[det_id]->SegmentEn[seg];
					if(En>0) 
					{
						MBTuple->FillColumn(column_id, En/EnergyUnit);
						WroteMBTuple=true;
					}
				}
			}

	// fill ntuple for CryostatEn
	if(AnyMBTupleCryostatEn)
		for(size_t cryo_id=0; cryo_id<cryostats.size(); cryo_id++)
		{
			G4int column_id = cryostats[cryo_id]->MBTupleCryostatEnId;
			if(column_id!=FIELD_NOT_SET)
			{
				size_t nb_of_dets=cryostats[cryo_id]->detectorIDs.size();
				G4double cryo_en=0;
				for(size_t det_nb=0; det_nb<nb_of_dets; det_nb++)
				{
					size_t det_id = cryostats[cryo_id]->detectorIDs[det_nb];
					if(detectors[det_id]->CoreEnAccu > 0)
						cryo_en += detectors[det_id]->CoreEnAccu;
				}
				if(cryo_en>0) 
				{
					MBTuple->FillColumn(column_id, cryo_en/EnergyUnit);
					WroteMBTuple=true;
				}
			}
		}

	// fill ntuple for main interaction
	if(MBTupleSegmentNbId!=FIELD_NOT_SET && Interaction)
	{
		MBTuple->FillColumn(MBTupleSegmentNbId,
				(G4double)InteractionSegmentNb);
		WroteMBTuple=true;
	}
	if(MBTupleSegmentCodeId!=FIELD_NOT_SET && Interaction)
	{
		G4int segmentCode
			= MiniBallSingleDetector::GetSegmentCode(InteractionSegmentNb);
		MBTuple->FillColumn(MBTupleSegmentCodeId, (G4double)segmentCode);
		WroteMBTuple=true;
	}
	if(MBTupleDetectorNbId!=FIELD_NOT_SET && Interaction)
	{
		MBTuple->FillColumn(MBTupleDetectorNbId,
				(G4double)InteractionDetectorNb);
		WroteMBTuple=true;
	}
	if(MBTupleCryostatNbId!=FIELD_NOT_SET && Interaction)
	{
		MBTuple->FillColumn(MBTupleCryostatNbId,
				(G4double)InteractionCryostatNb);
		WroteMBTuple=true;
	}
	if(MBTupleEventNbId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleEventNbId,(G4double)InteractionEventNb);
		WroteMBTuple=true;
	}
	if(MBTupleTraceNbId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleTraceNbId,(G4double)InteractionTraceNb);
		WroteMBTuple=true;
	}
	if(MBTupleRadiusId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleRadiusId, InteractionPositionInDet.rho());
		WroteMBTuple=true;
	}
	if(MBTuplePhiId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTuplePhiId, InteractionPositionInDet.phi());
		WroteMBTuple=true;
	}
	if(MBTupleZId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleZId, InteractionPositionInDet.z());
		WroteMBTuple=true;
	}
	if(MBTupleXId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleXId, InteractionPositionInDet.x());
		WroteMBTuple=true;
	}
	if(MBTupleYId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleYId, InteractionPositionInDet.y());
		WroteMBTuple=true;
	}
	if(MBTupleWorldXId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleWorldXId, InteractionPositionInWorld.x());
		WroteMBTuple=true;
	}
	if(MBTupleWorldYId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleWorldYId, InteractionPositionInWorld.y());
		WroteMBTuple=true;
	}
	if(MBTupleWorldZId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleWorldZId, InteractionPositionInWorld.z());
		WroteMBTuple=true;
	}
	if(MBTupleWorldThetaId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleWorldThetaId,
				InteractionPositionInWorld.theta());
		WroteMBTuple=true;
	}
	if(MBTupleWorldPhiId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleWorldPhiId,
				InteractionPositionInWorld.phi());
		WroteMBTuple=true;
	}
	if(MBTupleFIRadiusId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIRadiusId,
				FirstInteractionPositionInDet.rho());
		WroteMBTuple=true;
	}
	if(MBTupleFIPhiId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIPhiId, FirstInteractionPositionInDet.phi());
		WroteMBTuple=true;
	}
	if(MBTupleFIZId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIZId, FirstInteractionPositionInDet.z());
		WroteMBTuple=true;
	}
	if(MBTupleFIXId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIXId, FirstInteractionPositionInDet.x());
		WroteMBTuple=true;
	}
	if(MBTupleFIYId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIYId, FirstInteractionPositionInDet.y());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldXId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIWorldXId,
				FirstInteractionPositionInWorld.x());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldYId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIWorldYId,
				FirstInteractionPositionInWorld.y());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldZId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIWorldZId,
				FirstInteractionPositionInWorld.z());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldThetaId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIWorldThetaId,
				FirstInteractionPositionInWorld.theta());
		WroteMBTuple=true;
	}
	if(MBTupleFIWorldPhiId!=FIELD_NOT_SET && InteractionInEvent)
	{
		MBTuple->FillColumn(MBTupleFIWorldPhiId,
				FirstInteractionPositionInWorld.phi());
		WroteMBTuple=true;
	}

	if(MBTupleManagement==INTERN_MANAGEMENT) {
		if(WroteMBTuple){
			MBTuple->Fill();
		}
		else{
			MBTuple->Clear();
		}
	}
}

G4double MiniBallHistoManager::GetSegmentEn(size_t cryo_id,
		size_t det_nb, size_t seg_id)
{
	//converting from relative det_nb to absolut det_id
	if(cryo_id >= cryostats.size())
		return(0.);

	if(det_nb >= cryostats[cryo_id]->detectorIDs.size())
		return(0.);

	G4int det_id=cryostats[cryo_id]->detectorIDs[det_nb];

	return(GetSegmentEn(det_id, seg_id));
}

G4double MiniBallHistoManager::GetSegmentEn(size_t det_id, size_t seg_id)
{
	// check whether detector exists
	if(det_id >= detectors.size())
		return(0.);

	G4double value=0.;

	if(seg_id == CORE)
		value = detectors[det_id]->CoreEnAccu;
	else if(seg_id<(size_t)NB_OF_SEGMENTS)
		value = detectors[det_id]->SegmentEnAccu[seg_id];

	return(value);
}

G4double MiniBallHistoManager::GetAddbackEn(size_t cryo_id)
{
	//converting from relative det_nb to absolut det_id
	if(cryo_id >= cryostats.size())
		return(0.);

	size_t nb_of_dets=cryostats[cryo_id]->detectorIDs.size();
	G4double cryo_en=0;
	for(size_t det_id=0; det_id<nb_of_dets; det_id++)
		if(detectors[det_id]->CoreEnAccu > 0)
			cryo_en += detectors[det_id]->CoreEnAccu;

	return(cryo_en);
}



