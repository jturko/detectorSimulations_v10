#ifndef ISOTOPES
#define ISOTOPES

#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <iostream>
#include <vector>

#include "Element.hh"

int CalculateQ(Element**,double*,double*);

class Isotopes : Element {
public:
  Isotopes(const char* massFile = "./mass.dat", bool debug = false);
  ~Isotopes();
  Element* Search(char*);
  Element* Search(int, int);
  double Mass(int, int);
  double MassExcess(int, int);
  void ProjectileOneNeutronGain(Element*&, Element*&);
  void ProjectileOneNeutronLoss(Element*&, Element*&);
  void ProjectileTwoNeutronGain(Element*&, Element*&);
  void ProjectileTwoNeutronLoss(Element*&, Element*&);
  Element* Reaction(Element*, int, int);
  size_t NumberOfIsotopes() { return fIsotopeTable.size(); }

private:
  size_t GetIndex(int, int);
  std::vector<Element*> fIsotopeTable;
  static double fAmu;
  static bool fDebug;

  //ClassDef(Isotopes, 1);
};

#endif
