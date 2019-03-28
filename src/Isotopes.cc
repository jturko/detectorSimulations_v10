/****************************************************
  This file contains the code for handling the mass tables.
 ******************************************************/
#include "Isotopes.hh"
#include <fstream>

double Isotopes::fAmu = 931.494043;//error is 0.000080
bool Isotopes::fDebug = false;

//#define debug
//#define debug_all
//#define debug_transfer

Isotopes::Isotopes(const char* massFileName, bool debug) {
	fDebug = debug;

	char temp[256],*cptr;
	int i = 0;
	int n;
	int z;
	double emass;
	double dmass;
	std::ifstream massFile;
	massFile.open(massFileName,std::ios::in);

	while(!massFile.bad() && !massFile.eof()) {
		massFile>>z;
		massFile>>n;
		massFile>>emass;
		massFile>>dmass;
		massFile.width(6);
		massFile>>temp;
		cptr = temp;
		while( *cptr == ' ')	{
			cptr++;
		}
		while(isdigit(*cptr)) {
			cptr++;
		}
		fIsotopeTable.push_back(new Element(z,n,((double)(n+z)*fAmu)+emass/1000.,emass/1000.,dmass/1000.,cptr,i));
		++i;
		massFile.ignore(256,'\n');
		if(fDebug) {
			std::cout<<"read: "<<z<<" "<<n<<" "<<emass<<" "<<dmass<<" "<<temp<<"/"<<cptr<<std::endl;
		}
	}

	massFile.close();
}

Isotopes::~Isotopes() {
	for(size_t i = 0; i < fIsotopeTable.size(); ++i) {
		delete fIsotopeTable[i];
	}
}

size_t Isotopes::GetIndex(int z, int n) {
	size_t i = 0;
	for(i = 0; i < fIsotopeTable.size(); ++i) {
		if(fIsotopeTable[i]->Z() == z && fIsotopeTable[i]->N() == n) break;
		if(fDebug) {
			std::cout<<i<<". "<<fIsotopeTable[i]->Z()<<" "<<fIsotopeTable[i]->N()<<"("<<z<<", "<<n<<")"<<std::endl;
		}
	}

	return i;
}

double Isotopes::Mass(int z, int n) {
	if(fDebug) {
		std::cout<<"start Mass("<<z<<","<<n<<")"<<std::endl;
	}

	size_t i = GetIndex(z, n);
	if(i == NumberOfIsotopes()) {
		return -1.;
	}

	return fIsotopeTable[i]->Mass();
}

double Isotopes::MassExcess(int z, int n) {
	if(fDebug) {
		std::cout<<"start MassExcess("<<z<<","<<n<<")"<<std::endl;
	}

	size_t i = GetIndex(z, n);
	if(i == NumberOfIsotopes()) {
		return -1.;
	}

	return fIsotopeTable[i]->MassExcess();
}

Element* Isotopes::Search(char* c) {
	if(fDebug) {
		std::cout<<"start SearchTable("<<c<<")"<<std::endl;
	}
	int i;
	int n;
	char* cptr;
	char* tmp;
	double a;

	cptr = c;

	//check whether c starts with numbers
	if(isalpha(*cptr)) {
		std::cerr<<"Error in search_table: string "<<c<<" doesn't start with mass number"<<std::endl;
		exit(1);
	}

	//get number of digits in c and move cptr to first non-digit
	n = 0;

	while(isdigit(*cptr)) {
		n++;
		cptr++;
	}

	//copy letters from c to cptr and convert to lower (except for first letter)
	tmp = cptr;
	cptr = new char[strlen(cptr)+1];
	strcpy(cptr,tmp);
	i = 1;
	while (*(cptr+i)) {
		*(cptr+i)=tolower(*(cptr+i));
		i++;
	}

	//copy digits from c to tmp and convert to double
	tmp = new char[n+1];
	strncpy(tmp,c,n);
	tmp[n] = '\0';
	a = atof(tmp);

	i = 0;
	for(i = 0; i < (int) NumberOfIsotopes(); ++i) {
		if(strcmp(fIsotopeTable[i]->Name(),cptr)==0) {
			if(fDebug) {
				std::cout<<"Search::A(): "<<fIsotopeTable[i]->A()<<std::endl;
			}
			while(fIsotopeTable[i]->A() != a && i < (int) NumberOfIsotopes())	{
				if(fDebug) {
					std::cout<<a<<" != "<<fIsotopeTable[i]->A()<<" ("<<i<<")"<<std::endl;
				}
				++i;
			}
			break;
		}
	}

	delete[] cptr;
	delete[] tmp;

	if(i < (int) NumberOfIsotopes()) {
		return fIsotopeTable[i];
	} else {
		return NULL;
	}
}

Element* Isotopes::Search(int z, int n) {
	if(fDebug) {
		std::cout<<"start SearchTable("<<z<<", "<<n<<")"<<std::endl;
	}
	size_t i = GetIndex(z, n);
	if(i == NumberOfIsotopes()) {
		return NULL;
	}

	return fIsotopeTable[i];
}

void Isotopes::ProjectileOneNeutronGain(Element*& ejectile, Element*& recoil) {
	if(fDebug) {
		std::cout<<"called Isotopes::ProjectileOneNeutronGain("<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<")"<<std::endl;
	}
	ejectile = Reaction(ejectile,0,1);
	recoil   = Reaction(recoil,0,-1);
	if(fDebug) {
		std::cout<<"returned: "<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<std::endl;
	}
}

void Isotopes::ProjectileOneNeutronLoss(Element*& ejectile, Element*& recoil) {
	if(fDebug) {
		std::cout<<"called Isotopes::ProjectileOneNeutronLoss("<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<")"<<std::endl;
	}
	ejectile = Reaction(ejectile,0,-1);
	recoil   = Reaction(recoil,0,1);
	if(fDebug) {
		std::cout<<"returned: "<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<std::endl;
	}
}

void Isotopes::ProjectileTwoNeutronGain(Element*& ejectile, Element*& recoil) {
	if(fDebug) {
		std::cout<<"called Isotopes::ProjectileTwoNeutronGain("<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<")"<<std::endl;
	}
	ejectile = Reaction(ejectile,0,2);
	recoil   = Reaction(recoil,0,-2);
	if(fDebug) {
		std::cout<<"returned: "<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<std::endl;
	}
}

void Isotopes::ProjectileTwoNeutronLoss(Element*& ejectile, Element*& recoil) {
	if(fDebug) {
		std::cout<<"called Isotopes::ProjectileTwoNeutronLoss("<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<")"<<std::endl;
	}
	ejectile = Reaction(ejectile,0,-2);
	recoil   = Reaction(recoil,0,2);
	if(fDebug) {
		std::cout<<"returned: "<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<std::endl;
	}
}

Element* Isotopes::Reaction(Element* el, int protonChange, int neutronChange) {
	if(fDebug) {
		std::cout<<"called Isotopes::Reaction("<<el->A()<<el->Name()<<", "<<protonChange<<", "<<neutronChange<<")"<<std::endl;
	}
	int newProtonNumber, newNeutronNumber, index;

	if(protonChange == 0 && neutronChange == 0) {
		return el;
	}

	newProtonNumber = el->Z() + protonChange;
	newNeutronNumber = el->N() + neutronChange;
	index = el->Index();

	if(fDebug) {
		std::cout<<"newProtonNumber: "<<newProtonNumber<<", newNeutronNumber: "<<newNeutronNumber<<", index: "<<index<<std::endl;
	}

	if(protonChange >= 0 && neutronChange >= 0) {
		while(fIsotopeTable[index]->Z() != newProtonNumber || fIsotopeTable[index]->N() != newNeutronNumber) {
			index++;
			if(fDebug) {
				std::cout<<"increment index: "<<index<<std::endl;
			}
		}
	} else {
		while(fIsotopeTable[index]->Z() != newProtonNumber || fIsotopeTable[index]->N() != newNeutronNumber) {
			index--;
			if(fDebug) {
				std::cout<<"decrement index: "<<index<<std::endl;
			}
		}
	}

	if(fDebug) {
		std::cout<<"final index: "<<index<<" = "<<fIsotopeTable[index]->A()<<fIsotopeTable[index]->Name()<<std::endl;
	}

	return fIsotopeTable[index];
}

int CalculateQ(Element** particle,double* q,double* dq) {
	double dummy;
	/*********** calculate q-value *******************************/
	*q = -1.0*(particle[2]->MassExcess()+particle[3]->MassExcess());
	*q += (particle[0]->MassExcess()+particle[1]->MassExcess());
	//printf("CaluculateQ: q = %f + %f - %f - %f (%.0f%s + %.0f%s - %.0f%s - %.0f%s)\n",particle[0]->MassExcess(),particle[1]->MassExcess(),particle[2]->MassExcess(),particle[3]->MassExcess(),particle[0]->A(),particle[0]->Name(),particle[1]->A(),particle[1]->Name(),particle[2]->A(),particle[2]->Name(),particle[3]->A(),particle[3]->Name());
	dummy = particle[0]->DMass();
	*dq = dummy*dummy;
	dummy = particle[1]->DMass();
	*dq += dummy*dummy;
	dummy = particle[2]->DMass();
	*dq += dummy*dummy;
	dummy = particle[3]->DMass();
	*dq += dummy*dummy;
	*dq = (double)sqrt(*dq);

	return 0;
}

