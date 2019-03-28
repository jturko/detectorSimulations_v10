#ifndef MATERIAL
#define MATERIAL 1

#include <stdio.h>
#include <ctype.h>
#include <iostream>
#include <vector>

#include "Element.hh"

class Material {
	public:
		Material();
		Material(char*);
		Material(char*, bool);
		Material(int, char**, double*, bool);
		~Material();

		size_t NumberOfElements();
		bool Gaseous();
		char* Name();
		double Mass();
		Element* GetElement(size_t);
		double Fraction(size_t);
		double A(size_t);
		double Mass(size_t);
		int Z(size_t);
		char* Name(size_t);

	private:
		char* fCode;
		std::vector<Element*> fElements;
		std::vector<double> fFractions;
		bool fGaseous;
		double fMass;

		//ClassDef(Material, 1);
};

#endif
