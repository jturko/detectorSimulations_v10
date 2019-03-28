/************************************************************************
 * 
 *  handle histogramming
 * 
 ************************************************************************/

#include "MiniBallHistoGenerator.hh"

MiniBallHistogram::~MiniBallHistogram()
{
}

void MiniBallHistogram::Fill(G4double)
{
  std::cerr << "MiniBallHistogram: trying to use virtual class" << std::endl;
}

void MiniBallHistogram::SetBinContent(G4int, G4double)
{
  std::cerr << "MiniBallHistogram: trying to use virtual class" << std::endl;
}

void MiniBallHistogram::SetBinError(G4int, G4double)
{
  std::cerr << "MiniBallHistogram: trying to use virtual class" << std::endl;
}

G4double MiniBallHistogram::GetBinContent(G4int)
{
  std::cerr << "MiniBallHistogram: trying to use virtual class" << std::endl;
  return(-1.);
}

MiniBallTuple::~MiniBallTuple()
{
}

G4int MiniBallTuple::AddColumn(G4String, G4double)
{
  std::cerr << "MiniBallTuple: trying to use virtual class" << std::endl;
  return(-1);
}

void MiniBallTuple::FillColumn(G4int, G4double)
{
  std::cerr << "MiniBallTuple: trying to use virtual class" << std::endl;
}

G4int MiniBallTuple::AddLV(G4double)
{
  std::cerr << "MiniBallTuple: trying to use virtual class" << std::endl;
  return(-1);
}

void MiniBallTuple::FillLV(G4int, G4LorentzVector)
{
  std::cerr << "MiniBallTuple: trying to use virtual class" << std::endl;
}

G4int MiniBallTuple::AddArray(G4String, size_t, G4double)
{
  std::cerr << "MiniBallTuple: trying to use virtual class" << std::endl;
  return(-1);
}

void MiniBallTuple::FillArray(G4int, size_t, G4double)
{
  std::cerr << "MiniBallTuple: trying to use virtual class" << std::endl;
}

void MiniBallTuple::Fill()
{
  std::cerr << "MiniBallTuple: trying to use virtual class" << std::endl;
}

void MiniBallTuple::Clear()
{
  std::cerr << "MiniBallTuple: trying to use virtual class" << std::endl;
}

MiniBallHistoGenerator::~MiniBallHistoGenerator()
{
}

MiniBallHistogram* MiniBallHistoGenerator::NewHistogram1D(G4String,
					       G4int, G4double, G4double)
{
  std::cerr << "MiniBallHistoGenerator: trying to use virtual class" << std::endl;
  return(NULL);
}

MiniBallTuple* MiniBallHistoGenerator::NewTuple(G4String)
{
  std::cerr << "MiniBallHistoGenerator: trying to use virtual class" << std::endl;
  return(NULL);
}

void MiniBallHistoGenerator::Close()
{
  std::cerr << "MiniBallHistoGenerator: trying to use virtual class" << std::endl;
}
