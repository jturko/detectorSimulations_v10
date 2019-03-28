/************************************************************************
 * 
 *  handle histogramming using Root
 * 
************************************************************************/

#include "MiniBallRootGenerator.hh"

#include "TFile.h"

#include <string>

#include <stdio.h>

MiniBallRootHistogram::MiniBallRootHistogram(G4String histo_name, G4int bin,
					     G4double low, G4double high)
{
  hist = new TH1F(histo_name, histo_name, bin, low, high);
  hist_written=false;
}

MiniBallRootHistogram::~MiniBallRootHistogram()
{
  Close();
}

void MiniBallRootHistogram::Close()
{
  if(hist && !hist_written) {
    // write histogram to file
    hist->Write();
    hist_written=true;
  }
}

void MiniBallRootHistogram::Fill(G4double value)
{
  if(!hist)
    return;

  hist->Fill((Double_t)value);
}

void MiniBallRootHistogram::SetBinContent(G4int bin, G4double content)
{
  if(!hist)
    return;

  G4int nb_of_bins = hist->GetNbinsX();

  if(bin<0 || bin>nb_of_bins)
    return;

  hist->SetBinContent(bin+1, content);
}

void MiniBallRootHistogram::SetBinError(G4int bin, G4double error)
{
  if(!hist)
    return;

  G4int nb_of_bins = hist->GetNbinsX();

  if(bin<0 || bin>nb_of_bins)
    return;

  hist->SetBinError(bin+1, error);
}

G4double MiniBallRootHistogram::GetBinContent(G4int bin)
{
  if(!hist)
    return(0.);

  G4int nb_of_bins = hist->GetNbinsX();

  if(bin<0 || bin>nb_of_bins)
    return(0.);

  return(hist->GetBinContent(bin+1));
}

///////////////////////////////////////////////////////////////////////////////
// ROOT-tuple
///////////////////////////////////////////////////////////////////////////////

MiniBallRootTuple::MiniBallRootTuple(G4String tuple_name, G4int tree_size)
  : started_writing(false)
{
  rtree = new TTree(tuple_name, tuple_name);
  Long64_t GByte = 1000000000L;
  rtree->SetMaxTreeSize((Long64_t)tree_size*GByte);
  rtree->SetAutoSave();
  tree_written=false;
  tcaLV_init=false;
}

MiniBallRootTuple::~MiniBallRootTuple()
{
  Close();
}

void MiniBallRootTuple::Close()
{
  if(rtree && !tree_written) {
    rtree->GetCurrentFile()->cd();
    rtree->Write();
    tree_written=true;
  }
}

G4int MiniBallRootTuple::AddColumn(G4String field_name, G4double default_value)
{
  if(started_writing)
    {
      G4cerr << "MiniBallRootTuple: cannot create new column after "
	     << "started writing data" << std::endl;
      return(-1);
    }
  if(column_values.size()==ARRAY_IDS-1)
    {
      G4cerr << "MiniBallRootTuple: cannot create more than " << ARRAY_IDS
	     << " columns" << std::endl;
      return(-1);
    }
  G4cout << "adding field " << field_name << std::endl;  

  //Pointer for the new_entry
  Double_t* new_entry = new Double_t;
  //new_entry is filled with the default value
  *new_entry = default_value;
  //vector of pointers
  column_values.push_back(new_entry);
  //vector of default values for the method clear()
  column_default.push_back(default_value);
  //vector of strings for the field name
  column_names.push_back(field_name);

  //column_id of the new created column
  G4int column_id = column_values.size()-1;

  //new Branch in the ROOT-tree
  G4String field_type = field_name+"/D";
  rtree->Branch(field_name, column_values[column_id], field_type);

  return(column_id);
}

void MiniBallRootTuple::FillColumn(G4int column_id, G4double value)
{
  // check if column_id is valid
  if(column_id<0 || column_id>=(G4int)column_values.size())
    return;
  *(column_values[column_id])=value;
}


G4int MiniBallRootTuple::AddLV(G4double default_value)
{
  if(started_writing)
    {
      G4cerr << "MiniBallRootTuple: cannot create new column after "
	     << "started writing data" << std::endl;
      return(-1);
    }

  //check if TClonesArray already created
  if(!tcaLV_init) {
    // create TClonesArray
    tcaLV = new TClonesArray("TLorentzVector", 4);
    // add branch to tree
    rtree->Branch("tcaLV", "TClonesArray", &tcaLV, 32000, 0);
    tcaLV_init=true;
  }
  
  // check actual size of 'tcaLV'
  if(tcaLV->GetEntries() == tcaLV->GetSize())
    {
      G4cerr << "MiniBallRootTuple: cannot create more than " << tcaLV->GetSize()
	     << " TLorentzVectors -> increase size of 'tcaLV[" << tcaLV->GetSize() << "]'" << std::endl;
      return(-1);
    }

  // get entries in tcaLV
  Int_t new_tcaLV_index = tcaLV->GetEntries();

  // add new object
  new ( (*tcaLV)[new_tcaLV_index] ) TLorentzVector();

  // fill new object with default values
  ((TLorentzVector *)tcaLV->At(new_tcaLV_index))->SetPxPyPzE(default_value, default_value, default_value, default_value);

  G4cout << "added TLorentzVector 'tcaLV[" << new_tcaLV_index << "]'" << std::endl;  

  return(new_tcaLV_index);
}


void MiniBallRootTuple::FillLV(G4int G4LV_id, G4LorentzVector G4LV)
{
  // copy "G4LV" to TLorentzVector "TLV"
  TLorentzVector TLV;
  TLV.SetPxPyPzE(G4LV.px(), G4LV.py(), G4LV.pz(), G4LV.e());
  // set "tcaLV[G4LV_id]" to TLV 
  *( (TLorentzVector *)tcaLV->At(G4LV_id) ) = TLV;

}


G4int MiniBallRootTuple::AddArray(G4String array_name, size_t array_length,
				  G4double default_value)
{
  if(started_writing)
    {
      G4cerr << "MiniBallRootTuple: cannot create new column after "
	   << "started writing data" << std::endl;
      return(-1);
  }
  G4cout << "adding array " << array_name
	 << "[" << array_length << "]" << std::endl;  

  //Pointer for the new_entry
  Double_t* new_entry = new Double_t[array_length];
  //fill new_entry with default values
  for(size_t i=0; i<array_length; i++)
    new_entry[i] = default_value;
  //add pointer to vector
  array_values.push_back(new_entry);
  array_sizes.push_back(array_length);
  //vector of default values for the method clear()
  array_default.push_back(default_value);
  //vector of strings for the field name
  array_names.push_back(array_name);

  //array_id of the new created array
  G4int array_id = array_values.size()-1 + ARRAY_IDS;

  //new Branch in the ROOT-tree
  char array_type[500];
  sprintf(array_type, "%s[%d]/D", array_name.c_str(), (int)array_length);
  rtree->Branch(array_name, array_values[array_id-ARRAY_IDS], array_type);
  
  return(array_id);
}

void MiniBallRootTuple::FillArray(G4int array_id, size_t index,
				  G4double value)
{
  // check if array_id is valid
  if(array_id<ARRAY_IDS || array_id>=(G4int)array_values.size()+ARRAY_IDS)
    return;
  if(index>=array_sizes[array_id-ARRAY_IDS])
    return;
  Double_t* entry=array_values[array_id-ARRAY_IDS];
  entry[index]=value;
}

void MiniBallRootTuple::Clear()
{
  for(size_t column_id=0; column_id<column_values.size(); column_id++)
     *(column_values[column_id])=column_default[column_id];

  for(size_t array_nb=0; array_nb<array_values.size(); array_nb++)
    {
      Double_t* entry=array_values[array_nb];
      for(size_t index=0; index<array_sizes[array_nb]; index++)
        entry[index]=array_default[array_nb];
    }
}

void MiniBallRootTuple::Fill()
{
  rtree->Fill();
  Clear();
}

///////////////////////////////////////////////////////////////////////////////
// Manager for all kinds of ROOT-histogramming
///////////////////////////////////////////////////////////////////////////////

MiniBallRootGenerator::MiniBallRootGenerator(G4String root_file_name,
					     G4int max_tree_size_GByte)
  : RootFileName(root_file_name), max_tree_size(max_tree_size_GByte),
    root_file_opened(false), root_file_written(false), tf(NULL)
{;}

MiniBallRootGenerator::~MiniBallRootGenerator()
{
  Close();
}

MiniBallHistogram* MiniBallRootGenerator::NewHistogram1D(G4String histo_name,
		             G4int bin, G4double low, G4double high)
{
  if(!root_file_opened)
    {
      Open();
    }
  MiniBallRootHistogram* new_hist
    = new MiniBallRootHistogram(histo_name, bin, low, high);
  G4cout << "creating MiniBall Histogram " << histo_name << G4endl; 

  histograms.push_back(new_hist);

  return(new_hist);
}

MiniBallTuple* MiniBallRootGenerator::NewTuple(G4String tuple_name)
{
  if(!root_file_opened)
    {
      Open();
    }
  MiniBallRootTuple* new_tuple = new MiniBallRootTuple(tuple_name,
						       max_tree_size);

  G4cout << "creating MiniBall RootTree " << tuple_name << G4endl;

  tuples.push_back(new_tuple);

  return(new_tuple);
}

void MiniBallRootGenerator::Open()
{
  if(!root_file_opened)
    {
      tf = new TFile(RootFileName, "RECREATE");
      tf->cd();
      G4cout << "Opening the file "
	     << RootFileName
	     << " for collecting the data."
	     << G4endl;
      root_file_opened=true;
    }
}

void MiniBallRootGenerator::Close()
{
  if(!root_file_written)
      {
	for(size_t hist_id=0; hist_id<histograms.size(); hist_id++)
	  histograms[hist_id]->Close();

	for(size_t tuple_id=0; tuple_id<tuples.size(); tuple_id++)
	  tuples[tuple_id]->Close();

	if(tf)
	  {
	    tf->Close();
	    delete tf;
	  }
	root_file_written=true;
      }
}

