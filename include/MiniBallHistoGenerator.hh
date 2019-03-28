//////////////////////////////////////////////////////////////////////////////
/**  \file 
 *   \brief virtual handling of histogramming
 *
 *   \author Hans Boie
 */
////////////////////////////////////////////////////////////////////////////// 
#ifndef MiniBallHistoGenerator_h 
#define MiniBallHistoGenerator_h 

#include "globals.hh"
#include "G4LorentzVector.hh"

/// Virtual class for an histogramm
class MiniBallHistogram
{
public:
  virtual ~MiniBallHistogram();
  
public:
  /// Increment an histogramm at a given value.
  virtual void Fill(G4double x_value);

  /// set content of a given bin
  virtual void SetBinContent(G4int bin, G4double content);

  /// set content of a given bin
  virtual void SetBinError(G4int bin, G4double error);

  /// returns content of a given bin
  virtual G4double GetBinContent(G4int bin);
};

/// Virtual class for an ntuple.
class MiniBallTuple
{
public:
  virtual ~MiniBallTuple();

public:
  /// Add new column to ntuple with name field_name.
  /** Returns the identifier with which the column must be referred to.*/
  virtual G4int AddColumn(G4String field_name, G4double default_value=0);

  /// Fill a given value into a column of a ntuple. 
  /** The column must be referred to using the column identifier. */
  virtual void FillColumn(G4int column_id, G4double value);

  /// add and fill TLorentzVector
  virtual G4int AddLV(G4double default_value);
  virtual void FillLV(G4int G4LV_id, G4LorentzVector G4LV);
  
  /// Add new array to ntuple with name array_name and size array_size.
  /** Returns the identifier with which the array must be referred to.*/
    virtual G4int AddArray(G4String array_name, size_t array_size,
			   G4double default_value=0);

  /// Fill a given value into a array in a ntuple at a certain index number
  /** The array must be referred To using the array identifier.*/
  virtual void FillArray(G4int array_id, size_t index, G4double value);
  
  /// Fill all values into the ntuple.
  virtual void Fill();
  
  /// Clear all values, do not fill them.
  virtual void Clear();
};

/// Virtual class for the generation of histogramming
class MiniBallHistoGenerator 
{
public:
  virtual ~MiniBallHistoGenerator();

public:
  /// Generates a new 1-dim histogram. 
  /** MiniBallHistoGenerator will care for it,
      so the user should not delete the object!*/
  virtual MiniBallHistogram* NewHistogram1D(G4String histo_name,
		             G4int bin, G4double low, G4double high);

  // Generate a new ntuple.
  /** MiniBallHistoGenerator will care for it, 
      so the user should not delete the object!*/
  virtual MiniBallTuple* NewTuple(G4String tuple_name);

  /// Close file
  virtual void Close();
};
#endif
