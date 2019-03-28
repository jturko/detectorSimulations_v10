#include <stdio.h>
#include <ctype.h>
#include <iostream>

#include "Material.hh"

//#define debug

Material::Material() {
}

Material::Material(char* code_in) {
	int length = strlen(code_in);

	fCode = new char[length+1];

	for(int i=0;i<length;i++) {
		fCode[i] = toupper(code_in[i]);
	}
	fCode[length] = '\0';

	if(strstr(fCode,"PE") != 0) {
		if(isalpha(fCode[0])) {
			fGaseous = false;
			fElements.resize(2);
			fFractions.resize(2);

			fElements[0] = new Element((char*)"12C");
			fElements[1] = new Element((char*)"2H");
			fFractions[0] = 1./3.*fElements[0]->Mass()/(1./3.*fElements[0]->Mass()+2./3.*fElements[1]->Mass());
			fFractions[1] = 2./3.*fElements[1]->Mass()/(1./3.*fElements[0]->Mass()+2./3.*fElements[1]->Mass());
			//fMass = fElements[0]->Mass()+2.*fElements[1]->Mass();
			fMass = 2. * fElements[0]->Mass() + 4.*fElements[1]->Mass();
		} else {
			double ratio = atof(fCode);
			if(ratio < 0 || ratio > 1) {
				std::cout<<"Error, only ratios between 0 and 1 are allowed, tried to build "<<fCode<<" with deuteron/proton ratio = "<<ratio<<std::endl;
				exit(1);
			}
			if(ratio == 0) {
				fGaseous = false;
				fElements.resize(2);
				fFractions.resize(2);

				fElements[0] = new Element((char*)"12C");
				fElements[1] = new Element((char*)"1H");
				fFractions[0] = 1./3.*fElements[0]->Mass()/(1./3.*fElements[0]->Mass()+2./3.*fElements[1]->Mass());
				fFractions[1] = 2./3.*fElements[1]->Mass()/(1./3.*fElements[0]->Mass()+2./3.*fElements[1]->Mass());
				//fMass = fElements[0]->Mass()+2.*fElements[1]->Mass();
				fMass = 2. * fElements[0]->Mass() + 4.*fElements[1]->Mass();
			} else if(ratio == 1) {
				fGaseous = false;
				fElements.resize(2);
				fFractions.resize(2);

				fElements[0] = new Element((char*)"12C");
				fElements[1] = new Element((char*)"2H");
				fFractions[0] = 1./3.*fElements[0]->Mass()/(1./3.*fElements[0]->Mass()+2./3.*fElements[1]->Mass());
				fFractions[1] = 2./3.*fElements[1]->Mass()/(1./3.*fElements[0]->Mass()+2./3.*fElements[1]->Mass());
				//fMass = fElements[0]->Mass()+2.*fElements[1]->Mass();
				fMass = 2. * fElements[0]->Mass() + 4.*fElements[1]->Mass();
			} else {
				fGaseous = false;
				fElements.resize(3);
				fFractions.resize(3);

				fElements[0] = new Element((char*)"12C");
				fElements[1] = new Element((char*)"1H");
				fElements[2] = new Element((char*)"2H");
				fFractions[0] = 1./3.*fElements[0]->Mass()/(1./3.*fElements[0]->Mass()+(1-ratio)*2./3.*fElements[1]->Mass()+ratio*2./3.*fElements[2]->Mass());
				fFractions[1] = (1-ratio)*2./3.*fElements[1]->Mass()/(1./3.*fElements[0]->Mass()+(1-ratio)*2./3.*fElements[1]->Mass()+ratio*2./3.*fElements[2]->Mass());
				fFractions[2] = ratio*2./3.*fElements[2]->Mass()/(1./3.*fElements[0]->Mass()+(1-ratio)*2./3.*fElements[1]->Mass()+ratio*2./3.*fElements[2]->Mass());
				//fMass = fElements[0]->Mass()+(1-ratio)*2.*fElements[1]->Mass()+ratio*2.*fElements[2]->Mass();
				fMass = 2.* fElements[0]->Mass()+(1-ratio)*4.*fElements[1]->Mass()+ratio*4.*fElements[2]->Mass();
			}
		}
	} else if(strstr(fCode,"MY") != 0) {
		fGaseous = false;
		fElements.resize(3);
		fFractions.resize(3);

		fElements[0] = new Element((char*)"12C");
		fElements[1] = new Element((char*)"1H");
		fElements[2] = new Element((char*)"16O");
		fFractions[0] = 10./3.*fElements[0]->Mass()/(10./3.*fElements[0]->Mass()+8./3.*fElements[1]->Mass()+4./3.*fElements[2]->Mass());
		fFractions[1] =  8./3.*fElements[1]->Mass()/(10./3.*fElements[0]->Mass()+8./3.*fElements[1]->Mass()+4./3.*fElements[2]->Mass());
		fFractions[2] =  4./3.*fElements[2]->Mass()/(10./3.*fElements[0]->Mass()+8./3.*fElements[1]->Mass()+4./3.*fElements[2]->Mass());
		fMass = 10.*fElements[0]->Mass()+8.*fElements[1]->Mass()+4.*fElements[2]->Mass();
	} else {
		std::cerr<<"Material::Material(code): code = "<<fCode<<" ("<<code_in<<") not implemented yet!"<<std::endl;
	}
	std::cout << "creating Materials -----------------------" << std::endl;
#ifdef debug
	std::cout<<"Created instance of class Material("<<fCode<<") with "<<NumberOfElements()<<" elements:"<<std::endl;
	for(int i = 0; i < NumberOfElements(); ++i) {
		std::cout<<"    "<<i<<".: "<<A(i)<<Name(i)<<"/"<<fElements[i]->A()<<fElements[i]->Name()<<" (Z = "<<fElements[i]->Z()<<", N = "<<fElements[i]->N()<<")"<<std::endl;
	}
#endif
}

Material::Material(char* element, bool gas) {
	fGaseous = gas;
	fElements.resize(1);
	fFractions.resize(1);

	fElements[0] = new Element(element);
	fFractions[0] = 1;

	fMass = fElements[0]->Mass();

	fCode = new char[strlen(element)];
	strcpy(fCode,element);
}

Material::Material(int nofel, char* el[], double fr[], bool gas) {
	fGaseous = gas;
	fElements.resize(nofel);
	fFractions.resize(nofel);

	double sumOfFractions = 0.;
	fMass = 0;

	for(int i = 0; i < nofel; ++i) {
		fElements[i] = new Element(el[i]);
		fFractions[i] = fr[i];
		sumOfFractions += fr[i];
	}

	for(int i = 0; i < nofel; ++i) {
		fFractions[i] /= sumOfFractions;
		fMass += fElements[i]->Mass()*fFractions[i];
	}

	fCode = new char[3];
	fCode[0] = 'X';
	fCode[1] = 'X';
	fCode[2] = '\0';
}

Material::~Material() {
	for(size_t i = 0; i < NumberOfElements(); ++i) {
		if(fElements[i] != NULL) {
			delete fElements[i];
		}
	}

	if(fCode != NULL) {
		delete[] fCode;
	}
}

size_t Material::NumberOfElements() {
	return fElements.size();
}

bool Material::Gaseous() {
	return fGaseous;
}

char* Material::Name() {
	return fCode;
}

double Material::Mass() {
	return fMass;
}

Element* Material::GetElement(size_t el) {
	if(el < NumberOfElements()) return fElements[el];
	else return NULL;
}

double Material::Fraction(size_t el) {
	if(el < NumberOfElements()) return fFractions[el];
	return 0.;
}

double Material::A(size_t el) {
	if(el < NumberOfElements()) return fElements[el]->A();
	return 0.;
}

double Material::Mass(size_t el) {
	if(el < NumberOfElements()) return fElements[el]->Mass();
	return 0.;
}

int Material::Z(size_t el) {
	if(el < NumberOfElements()) return fElements[el]->Z();
	return -1;
}

char* Material::Name(size_t el) {
	if(el < NumberOfElements()) return fElements[el]->Name();
	return NULL;
}
