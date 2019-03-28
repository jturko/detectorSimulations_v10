#ifndef ELEMENT
#define ELEMENT

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <iostream>
#include <string>
#include <cstring>

class Element {
public:
  Element();
  Element(char*);
  Element(int, int, double, double, double, char*, int);
  ~Element();
  
  int Z();
  int N();
  double A();
  double Mass();
  double MassExcess();
  double DMass();
  char* Name();
  int Index();
  
private:
  int GetCharge(char*);
  double GetMass(int);
  
  double fMass;
  int fProtonNumber;
  int fNeutronNumber;
  double fMassExcess;
  double fDeltaMass;
  char* fCode;
  int fIndex;
  bool fNatural;

  //ClassDef(Element, 1);
};

#endif
