///////////////////////////////////////////////////////////////////////////////
/** \file
 *  \brief handle histogramming using ROOT
 *
 *  \author Hans Boie
 */
///////////////////////////////////////////////////////////////////////////////

#ifndef MiniBallRootGenerator_h
#define MiniBallRootGenerator_h 

#include "MiniBallHistoGenerator.hh"

#include <fstream>
#include <vector>

#include "TTree.h"
#include "TH1.h"
#include "G4LorentzVector.hh"
#include "TLorentzVector.h"
#include "TClonesArray.h"

// tree size in gigabyte
#define DEFAULT_TREE_SIZE 5

/// ROOT-histogramm. 
class MiniBallRootHistogram : public MiniBallHistogram
{
public:
  MiniBallRootHistogram(G4String hist_name, G4int nb_of_bins,
			G4double low, G4double high);
  ~MiniBallRootHistogram();
  
public:
  /// Increment an histogramm at a given x value.
  void Fill(G4double x_value);

  /// sets content of a given bin 
  /// bin numbers may range from 0 to nb_of_bins-1
  void SetBinContent(G4int bin, G4double content);

  /// sets error of a given bin 
  /// bin numbers may range from 0 to nb_of_bins-1
  void SetBinError(G4int bin, G4double error);

  /// returns content of a given bin 
  /// bin numbers may range from 0 to nb_of_bins-1
  G4double GetBinContent(G4int bin);

  /// Writes the histogramm to the file
  void Close();

private:
  TH1F* hist;

  G4bool hist_written;
};

#define TUPLE_FILE_EXT ".tup.dat"
#define ARRAY_IDS      10000

/// ROOT-Tuple
class MiniBallRootTuple : public MiniBallTuple
{
public:
  /// generate new MiniBallRootTuple
  MiniBallRootTuple(G4String tuple_name,
		    G4int tree_size_GByte=DEFAULT_TREE_SIZE);
  /**< the sieze of the tree is given in values of gigabyte */
  ~MiniBallRootTuple();
  
public:
  /// Add new column to ntuple with name field_name.
  G4int AddColumn(G4String field_name, G4double default_value=0);
  /**<  Returns the identifier with which the column must be referred to.*/
  
  /// Fill a given value into a column of a ntuple. 
  void FillColumn(G4int column_id, G4double value);
  /**< The column must be referred to using the column identifier.*/

  /// add and fill TLorentzVector 
  G4int AddLV(G4double default_value);
  void FillLV(G4int G4LV_id, G4LorentzVector G4LV);

  /// Add new array to ntuple with name array_name and size array_size.
  /** Returns the identifier with which the array must be referred to.*/
  G4int AddArray(G4String array_name, size_t array_size,
		 G4double default_value=0);

  /// Fill a given value into a array in a ntuple at a certain index number
  /** The array must be referred To using the array identifier.*/
  void FillArray(G4int array_id, size_t index, G4double value);

  /// Fill all values into the ntuple.
  void Fill();

  /// Clear all values, do not fill them.
  void Clear();

  /// Write the Tuple to the file
  void Close();

private:
  TTree* rtree;
  
  std::vector<Double_t*> column_values;
  std::vector<Double_t> column_default;
  std::vector<G4String> column_names;

  TClonesArray* tcaLV;
  G4bool tcaLV_init;

  std::vector<Double_t*> array_values;
  std::vector<size_t> array_sizes;
  std::vector<Double_t> array_default;
  std::vector<G4String> array_names;

  G4bool started_writing;
  G4bool tree_written;
};

/// Generator for all kinds of  ROOT histogramming
class MiniBallRootGenerator : public MiniBallHistoGenerator 
{
public:
  MiniBallRootGenerator(G4String root_file_name,
			G4int max_tree_size_GByte=DEFAULT_TREE_SIZE);
  ~MiniBallRootGenerator();

public:
  /// Generate a new 1-dim histogram
  MiniBallHistogram* NewHistogram1D(G4String histo_name,
				    G4int bin, G4double low, G4double high);
  
  /// Generate a new ntuple
  MiniBallTuple* NewTuple(G4String tuple_name);

  /// Close file
  void Close();
  /// Opens a file
  void Open();
  
private:
  G4String RootFileName;
  G4int max_tree_size;

  G4bool root_file_opened;
  G4bool root_file_written;

  TFile* tf;

  std::vector<MiniBallRootHistogram*> histograms;
  std::vector<MiniBallRootTuple*> tuples;
};
#endif




