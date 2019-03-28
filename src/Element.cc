#include "Element.hh"

//#define debug

Element::Element() {
}

Element::Element(char* el) {
	int length = strlen(el);
	if(length == 0) {
		std::cerr<<"ERROR, tried to create Element with empty string!"<<std::endl;
		exit(1);
	}

	int first_digit = length;
	int first_letter = length;
	int code_length = 0;
	int i;

	if(isdigit(el[0])) {
		first_digit = 0;
	} else if(isalpha(el[0])) {
		first_letter = 0;
	} else {
		std::cerr<<"ERROR, tried to create Element with string that doesn't start with an alphanumerial ("<<el<<")!"<<std::endl;
		exit(1);
	}

#ifdef debug
	std::cout<<"first_digit: "<<first_digit<<" and first_letter: "<<first_letter<<std::endl;
#endif

	for(i = 0; i < length; ++i) {
		if(isdigit(el[i]) && first_letter == 0) { //found first digit in a string that started with letters
			first_digit = i;
			code_length = i + 1;
			break;
		}
		if(isalpha(el[i]) && first_digit == 0) { //found first letter in a string that started with digits
			first_letter = i;
			code_length = length - i + 1;
			break;
		}
	}

	if(i == length) {
		if(first_digit == 0) { //string contains only digits
			std::cerr<<"ERROR, tried to create Element with string that contains only digits ("<<el<<")!"<<std::endl;
			exit(1);
		}
		if(first_letter == 0) { //string contains only letters
			fCode = new char[length+1];
			strcpy(fCode, el);

			fProtonNumber = GetCharge(fCode);
			fMass = GetMass(fProtonNumber);
			fNatural = true;
		}
	} else {
		fCode = new char[code_length];
		strncpy(fCode, el+first_letter, code_length - 1);
		fCode[code_length - 1] = '\0';

		fProtonNumber = GetCharge(fCode);
		fNeutronNumber = atoi(el+first_digit) - fProtonNumber;
		fMassExcess = 0.;
		fMass = atof(el+first_digit)*931.494034+fMassExcess;
		fNatural = false;
	}

#ifdef debug
	std::cout<<"Created instance of class element("<<el<<"): "<<flush<<A()<<Name()<<std::endl;
	std::cout<<"mass = "<<fMass<<", Z = "<<fProtonNumber<<", N = "<<fNeutronFumber<<", mass_excess = "<<fMassExcess<<", index = "<<fIndex<<std::endl;
#endif
}

Element::Element(int z,int n,double m,double m_ex,double dm,char* cd,int ind) {
	fProtonNumber = z;
	fNeutronNumber = n;
	fMass = m;
	fMassExcess = m_ex;
	fDeltaMass = dm;
	fCode = new char[strlen(cd)+1];
	strcpy(fCode,cd);
	fIndex = ind;
	fNatural = false;
#ifdef debug
	std::cout<<"Created instance of class element("<<z<<", "<<n<<", "<<m<<", "<<m_ex<<", "<<dm<<", "<<cd<<", "<<ind<<"): "<<flush<<A()<<Name()<<std::endl;
	std::cout<<"mass = "<<fMass<<", Z = "<<fProtonNumber<<", N = "<<fNeutronNumber<<", mass_excess = "<<fMassExcess<<", index = "<<fIndex<<std::endl;
#endif
}

Element::~Element() {
	if(fCode != NULL) {
		delete[] fCode;
	}
}

int Element::GetCharge(char* code) {
	char codes[105][3] = {"H","HE","LI","BE","B","C","N","O","F","NE","NA","MG","AL","SI","P","S","CL","AR","K","CA","SC","TI","V","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR","RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","F","XE","CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU","HF","TA","W","RE","OS","IR","PT","AU","HG","TI","PB","BI","PO","AT","RN","FR","RA","AC","TH","PA","U","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR","RF","HA"};

	int length = strlen(code);

	char* tmp = new char[length+1];

	for(int i = 0; i < length; ++i) {
		tmp[i] = toupper(code[i]);
	}
	tmp[length] = '\0';

	for(int i=0;i<105;i++) {
		if(strcmp(tmp,codes[i]) == 0) {
			delete[] tmp;
			return i+1;
		}
	}

	delete[] tmp;
	return 0;
}

double Element::GetMass(int proton_number) {
	double mass[106] = {1,1.00797,4.0026,6.939,9.0122,10.818,12.01115,14.0067,15.9994,18.99984,20.183,22.9898,24.312,26.9815,28.086,30.9738,32.064,35.453,39.948,39.102,40.08,44.956,47.90,50.942,51.996,54.938,55.847,58.933,58.71,63.54,65.37,69.72,72.59,74.922,78.96,79.909,83.80,85.47,87.62,88.909,91.22,92.906,95.94,98.,101.07,102.905,106.4,107.87,112.4,114.82,118.69,121.75,127.60,126.904,131.3,132.905,137.34,138.91,140.12,140.907,144.24,146.,150.35,151.96,157.25,158.924,162.50,164.93,167.26,168.934,173.04,174.97,178.49,180.948,183.85,186.2,190.2,192.2,195.09,196.967,200.59,204.37,207.19,208.98,210.,210.,222.,223.,226.,227.,232.038,231.,238.03,237.,242.,243.,247.,247.,249.,254.,253.,256.,254.,259.,257.,260.};

	if(0 <= proton_number && proton_number < 106) {
		return mass[proton_number];
	} else {
		return -1.;
	}
}

int Element::Z() {
	return fProtonNumber;
}

int Element::N() {
	return fNeutronNumber;
}

double Element::A() {
	if(fNatural) {
#ifdef debug
		std::cout<<"Element::A(): natural"<<std::endl;
#endif
		return fMass;
	} else {
#ifdef debug
		std::cout<<"Element::A(): not natural"<<std::endl;
#endif
		return ((double) (fProtonNumber+fNeutronNumber));
	}
}

double Element::Mass() {
	if(fNatural) {
		return fMass*931.494034;
	} else {
		return fMass;
	}
}

double Element::MassExcess() {
	return fMassExcess;
}

double Element::DMass() {
	return fDeltaMass;
}

char* Element::Name() {
	return fCode;
}

int Element::Index() {
	return fIndex;
}

