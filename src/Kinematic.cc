#include "Kinematic.hh"
#include <fstream>

//#define debug
//#define debug_orbits
//#define debug_irma
//#define debug_results
//#define e_strag_debug
//#define debug_orbits_results
//#define debug_vsx
//#define debug_VsX
//#define debug_VsE
//#define cmangle_debug

double lab2f(double ex, double cm, double Emean, Element** particle) {
	double p2[3],p3[3];
	double theta;
	double cm_unit[3];
	double e2=0;
	double lab;
	theta = cm * TMath::Pi() / 180.0;
	cm_unit[2] = cos(theta);
	cm_unit[1] = 0.0;
	cm_unit[0] = sin(theta);
	scatter(p2,p3,Emean,e2,ex,cm_unit,particle);
	lab=acos(p2[2]/norm(p2))*180.0/TMath::Pi();
#ifdef debug_orbits_results
	std::cout<<"lab2f("<<ex<<", "<<cm<<", "<<Emean<<", ...): "<<lab<<" = acos("<<p2[2]<<"/"<<norm(p2)<<")*180.0/TMath::Pi()"<<std::endl;
#endif
	return lab;
}

double lab3f(double ex, double cm, double Emean, Element** particle) {
	double p2[3],p3[3];
	double theta;
	double cm_unit[3];
	double e2=0;
	double lab;
		theta = cm * TMath::Pi() / 180.0;
	cm_unit[2] = cos(theta);
	cm_unit[1] = 0.0;
	cm_unit[0] = sin(theta);
	scatter(p2,p3,Emean,e2,ex,cm_unit,particle);
	lab=acos(p3[2]/norm(p3))*180.0/TMath::Pi();
#ifdef debug_orbits_results
	std::cout<<"lab2f("<<ex<<", "<<cm<<", "<<Emean<<", ...): "<<lab<<" = acos("<<p3[2]<<"/"<<norm(p3)<<")*180.0/TMath::Pi()"<<std::endl;
#endif
	return lab;
}

double en2f (double ex, double cm, double Emean, Element** particle) {
	double p2[3],p3[3];
	double theta;
	double cm_unit[3];
	double e2=0;
	double en;
	theta = cm * TMath::Pi() / 180.0;
	cm_unit[2] = cos(theta);
	cm_unit[1] = 0.0;
	cm_unit[0] = sin(theta);
	scatter(p2,p3,Emean,e2,ex,cm_unit,particle);
	en=energy(norm(p2),particle[2]->Mass()+e2);
	return en;
}

double en3f (double ex, double cm, double Emean, Element** particle) {
	double p2[3],p3[3];
	double theta;
	double cm_unit[3];
	double e2=0;
	double en;
	theta = cm * TMath::Pi() / 180.0;
	cm_unit[2] = cos(theta);
	cm_unit[1] = 0.0;
	cm_unit[0] = sin(theta);
	scatter(p2,p3,Emean,e2,ex,cm_unit,particle);
	en=energy(norm(p3),particle[3]->Mass()+ex);
	return en;
}


// ------------------------------------------------------------------------------

int Kinematic::orbits(Element* projectile, Element* target, Element* recoil, Element* ejectile, double energy, double cm_start, double cm_stop, double cm_step, double ex_start, double ex_stop, double ex_step, int code, double**& results_en, double**& results_ang) {
	double q,dq,e_cm,Ethr;
	double cm, ex;
	double lab2,lab3,en2,en3;

	Element* particle[4] = {target, projectile, ejectile, recoil};

	Ethr=0;
	ex=0;
	//load_table();  // load mass table

	CalculateQ(particle,&q,&dq);

	//
	e_cm = EnergyCm(energy,particle) - projectile->Mass() - target->Mass();


	// check for threshold
	if (e_cm < -q)  {
		printf("# Energy below threshold\n");
		exit(1);
	} 
	else if (q < -1.e-10)  {
		// find threshold
		double hEtemp,lEtemp,e_cmt;
		hEtemp=energy;
		lEtemp=0;
		do {
			Ethr=(hEtemp+lEtemp)/2;
			e_cmt = EnergyCm(Ethr,particle) - projectile->Mass() - target->Mass();
			if (e_cmt < -q) {
				lEtemp=Ethr;
			}
			else {
				hEtemp=Ethr;
			}
		} while ( fabs(e_cmt + q)/(-q) > 1e-5);
#ifdef debug_orbits
		printf("# Threshold = %f,  thr/A = %f \n",Ethr,Ethr/(projectile->Z() + projectile->N()));
#endif
	}

#ifdef debug_orbits
	printf("# Beam Energy = %f MeV\n",energy);
#endif

	// output

	results_en = new double*[((int)((ex_stop - ex_start)/ex_step))+1];
	results_ang = new double*[((int)((ex_stop - ex_start)/ex_step))+1];

	for(int i=0;i<((int)(ex_stop - ex_start)/ex_step)+1;i++) {
		results_en[i] = new double[((int)((cm_stop - cm_start)/cm_step))+1];
		results_ang[i] = new double[((int)((cm_stop - cm_start)/cm_step))+1];
	}

#ifdef debug_orbits
	std::cout<<"created new arrays results_en and results_ang with dimensions ["<<((int)((ex_stop - ex_start)/ex_step))+1<<"]["<<((int)((cm_stop - cm_start)/cm_step))+1<<"]"<<std::endl;
#endif

	for(int i=0;i<((int)(ex_stop - ex_start)/ex_step)+1;i++) {
		for(int j=0;j<((int)(cm_stop - cm_start)/cm_step)+1;j++) {
			results_en[i][j] = -1;
			results_ang[i][j] = -1;
		}
	}


	if(e_cm+q < ex_stop)  {
		ex_stop=e_cm+q;
		printf("# Readjusted stop value to %f, since below threshold.\n",ex_stop);
	}

	if(code == 0) {
#ifdef debug_orbits
		std::cout<<"#code 0: plot of lab energy as function of lab angles for ejectile"<<std::endl;
#endif
		for(ex=ex_start;ex<=ex_stop;ex +=ex_step) {
			for(cm=cm_start;cm<=cm_stop;cm +=cm_step) {
				lab2=lab2f(ex, cm, energy, particle);
				en2 = en2f(ex, cm, energy, particle);
				if(!(lab2>0) && !(lab2<0)) {
					lab2 = -1;
				}
#ifdef debug_orbits_results
				std::cout<<"lab2f("<<ex<<", "<<cm<<", "<<energy<<"): "<<lab2<<" "<<en2<<std::endl;
#endif
#ifdef debug_orbits
				std::cout<<ex<<"/"<<cm<<": results_en["<<(int)((ex - ex_start)/ex_step)<<"]["<<(int)((cm - cm_start)/cm_step)<<"] = "<<flush;
#endif
				results_en[(int)((ex - ex_start)/ex_step)][(int)((cm - cm_start)/cm_step)] = en2;
				results_ang[(int)((ex - ex_start)/ex_step)][(int)((cm - cm_start)/cm_step)] = lab2/180.*TMath::Pi();
#ifdef debug_orbits
				std::cout<<en2<<", results_ang["<<(int)((ex - ex_start)/ex_step)<<"]["<<(int)((cm - cm_start)/cm_step)<<"] = "<<lab2/180.*TMath::Pi()<<std::endl;
#endif
			}
		}
	}
	else if(code == 1) {
#ifdef debug_orbits
		std::cout<<"#code 1: plot of lab energy as function of lab angles for recoil"<<std::endl;
#endif
		for(ex=ex_start;ex<=ex_stop;ex +=ex_step) {
			for(cm=cm_start;cm<=cm_stop;cm +=cm_step) {
				lab3=lab3f(ex, cm, energy, particle);
				en3 = en3f(ex, cm, energy, particle);
				if(!(lab3>0) && !(lab3<0)) {
					lab3 = -1;
				}
#ifdef debug_orbits_results
				std::cout<<"lab3f("<<ex<<", "<<cm<<", "<<energy<<"): "<<lab3<<" "<<en3<<std::endl;
#endif
#ifdef debug_orbits
				std::cout<<ex<<"/"<<cm<<": results_en["<<(int)((ex - ex_start)/ex_step)<<"]["<<(int)((cm - cm_start)/cm_step)<<"] = "<<flush;
#endif
				results_en[(int)((ex - ex_start)/ex_step)][(int)((cm - cm_start)/cm_step)] = en3;
				results_ang[(int)((ex - ex_start)/ex_step)][(int)((cm - cm_start)/cm_step)] = lab3/180.*TMath::Pi();
#ifdef debug_orbits
				std::cout<<en3<<", results_ang["<<(int)((ex - ex_start)/ex_step)<<"]["<<(int)((cm - cm_start)/cm_step)<<"] = "<<lab3/180.*TMath::Pi()<<std::endl;
#endif
			}
		}
	}
	else if(code == 2) {
		std::cout<<"#code 2: plot of lab angles as function of cm angles"<<std::endl;
		for( ex = ex_start; ex <= ex_stop; ex += ex_step ) {
			for( cm = cm_start; cm <= cm_stop; cm += cm_step) {
				lab2=lab2f(ex, cm, energy, particle); 
				printf("%6.2f %12.5e \n", cm, lab2);
			}
		}
		for( ex = ex_start; ex <= ex_stop; ex += ex_step ) {
			for( cm = cm_start; cm <= cm_stop; cm += cm_step) {
				lab3=lab3f(ex, cm, energy, particle); 
				printf("%6.2f %12.5e \n", cm, lab3);
			}
		}
	}
	else if(code == 3) {
		// plot of lab energies as function of cm angles
		for( ex = ex_start; ex <= ex_stop; ex += ex_step ) {
			for( cm = cm_start; cm <= cm_stop; cm += cm_step) {
				en2=en2f(ex, cm, energy, particle); 
				printf("%6.2f %12.5e \n", cm, en2);
			}
		}
		for( ex = ex_start; ex <= ex_stop; ex += ex_step ) {
			for( cm = cm_start; cm <= cm_stop; cm += cm_step) {
				en3=en3f(ex, cm, energy, particle); 
				printf("%6.2f %12.5e \n", cm, en3);
			}
		}
	}

	return 0;
}

int Kinematic::orbits(Element* projectile, Element* target, Element* recoil, Element* ejectile, double energy, double cm_start, double cm_stop, double cm_step, double ex, int code, double*& results_en, double*& results_ang) {
	double q,dq,e_cm,Ethr;
	double cm;
	double lab2,lab3,en2,en3;

	Element* particle[4] = {target, projectile, ejectile, recoil};

	Ethr=0;
	//load_table();  // load mass table

	CalculateQ(particle,&q,&dq);

	//
	e_cm = EnergyCm(energy,particle) - projectile->Mass() - target->Mass();

	// check for threshold
	if (e_cm < -q)  {
		printf("# Energy below threshold\n");
		exit(1);
	} 
	else if (q < -1.e-10)  {
		// find threshold
		double hEtemp,lEtemp,e_cmt;
		hEtemp=energy;
		lEtemp=0;
		do {
			Ethr=(hEtemp+lEtemp)/2;
			e_cmt = EnergyCm(Ethr,particle) - projectile->Mass() - target->Mass();
			if (e_cmt < -q) {
				lEtemp=Ethr;
			}
			else {
				hEtemp=Ethr;
			}
		} while ( fabs(e_cmt + q)/(-q) > 1e-5);
		printf("# Threshold = %f,  thr/A = %f \n",Ethr,Ethr/(projectile->Z() + projectile->N()));
	}

#ifdef debug_orbits
	printf("# Beam Energy = %f MeV\n",energy);
#endif

	// output

	results_en = new double[((int)((cm_stop - cm_start)/cm_step))+1];
	results_ang = new double[((int)((cm_stop - cm_start)/cm_step))+1];

#ifdef debug_orbits
	std::cout<<"created new arrays results_en and results_ang with dimensions ["<<((int)((cm_stop - cm_start)/cm_step))+1<<"]"<<std::endl;
#endif

	for(int j=0;j<((int)(cm_stop - cm_start)/cm_step)+1;j++) {
		results_en[j] = -1;
		results_ang[j] = -1;
	}


	if(e_cm+q < ex)  {
		ex=e_cm+q;
		printf("# Readjusted excitation energy value to %f, since below threshold.\n",ex);
	}

	if(code == 0) {
#ifdef debug_orbits
		std::cout<<"#code 0: plot of lab energy as function of lab angles for ejectile"<<std::endl;
#endif
		for(cm=cm_start;cm<=cm_stop;cm +=cm_step) {
			lab2=lab2f(ex, cm, energy, particle);
			en2 = en2f(ex, cm, energy, particle);
			if(!(lab2>0) && !(lab2<0)) {
				lab2 = -1;
			}
#ifdef debug_orbits_results
			std::cout<<"lab2f("<<ex<<", "<<cm<<", "<<energy<<"): "<<lab2<<" "<<en2<<std::endl;
#endif
#ifdef debug_orbits
			std::cout<<ex<<"/"<<cm<<": results_en["<<(int)((cm - cm_start)/cm_step)<<"] = "<<flush;
#endif
			results_en[(int)((cm - cm_start)/cm_step)] = en2;
			results_ang[(int)((cm - cm_start)/cm_step)] = lab2/180.*TMath::Pi();
#ifdef debug_orbits
			std::cout<<en2<<", results_ang["<<(int)((cm - cm_start)/cm_step)<<"] = "<<lab2/180.*TMath::Pi()<<std::endl;
#endif
		}
	}
	else if(code == 1) {
#ifdef debug_orbits
		std::cout<<"#code 1: plot of lab energy as function of lab angles for recoil"<<std::endl;
#endif
		for(cm=cm_start;cm<=cm_stop;cm +=cm_step) {
			lab3=lab3f(ex, cm, energy, particle);
			en3 = en3f(ex, cm, energy, particle);
			if(!(lab3>0) && !(lab3<0)) {
				lab3 = -1;
			}
#ifdef debug_orbits_results
			std::cout<<"lab3f("<<ex<<", "<<cm<<", "<<energy<<"): "<<lab3<<" "<<en3<<std::endl;
#endif
#ifdef debug_orbits
			std::cout<<ex<<"/"<<cm<<": results_en["<<(int)((cm - cm_start)/cm_step)<<"] = "<<flush;
#endif
			results_en[(int)((cm - cm_start)/cm_step)] = en3;
			results_ang[(int)((cm - cm_start)/cm_step)] = lab3/180.*TMath::Pi();
#ifdef debug_orbits
			std::cout<<en3<<", results_ang["<<(int)((cm - cm_start)/cm_step)<<"] = "<<lab3/180.*TMath::Pi()<<std::endl;
#endif
		}
	}
	else if(code == 2) {
		std::cout<<"#code 2: plot of lab angles as function of cm angles"<<std::endl;
		for( cm = cm_start; cm <= cm_stop; cm += cm_step)  {
			lab2=lab2f(ex, cm, energy, particle); 
			printf("%6.2f %12.5e \n", cm, lab2);
		}
		for( cm = cm_start; cm <= cm_stop; cm += cm_step)  {
			lab3=lab3f(ex, cm, energy, particle); 
			printf("%6.2f %12.5e \n", cm, lab3);
		}
	}
	else if(code == 3) {
		// plot of lab energies as function of cm angles
		for( cm = cm_start; cm <= cm_stop; cm += cm_step)  {
			en2=en2f(ex, cm, energy, particle); 
			printf("%6.2f %12.5e \n", cm, en2);
		}
		for( cm = cm_start; cm <= cm_stop; cm += cm_step)  {
			en3=en3f(ex, cm, energy, particle); 
			printf("%6.2f %12.5e \n", cm, en3);
		}
	}

	return 0;
}

int Kinematic::orbits(Element* projectile, Element* target, Element* recoil, Element* ejectile, double energy, double cm_start, double cm_stop, double cm_step, double ex, int code, std::vector<double>& results_en, std::vector<double>& results_ang) {
	double q,dq,e_cm,Ethr;
	double cm;
	double lab2,lab3,en2,en3;

	Element* particle[4] = {target, projectile, ejectile, recoil};

	Ethr=0;
	//load_table();  // load mass table

	CalculateQ(particle,&q,&dq);

	//
	e_cm = EnergyCm(energy,particle) - projectile->Mass() - target->Mass();

	// check for threshold
	if (e_cm < -q)  {
		printf("# Energy below threshold\n");
		exit(1);
	} 
	else if (q < -1.e-10)  {
		// find threshold
		double hEtemp,lEtemp,e_cmt;
		hEtemp=energy;
		lEtemp=0;
		do {
			Ethr=(hEtemp+lEtemp)/2;
			e_cmt = EnergyCm(Ethr,particle) - projectile->Mass() - target->Mass();
			if (e_cmt < -q) {
				lEtemp=Ethr;
			}
			else {
				hEtemp=Ethr;
			}
		} while ( fabs(e_cmt + q)/(-q) > 1e-5);
		printf("# Threshold = %f,  thr/A = %f \n",Ethr,Ethr/(projectile->Z() + projectile->N()));
	}

#ifdef debug_orbits
	printf("# Beam Energy = %f MeV\n",energy);
#endif

	// output

	results_en.resize(((int)((cm_stop - cm_start)/cm_step))+1);
	results_ang.resize(((int)((cm_stop - cm_start)/cm_step))+1);

#ifdef debug_orbits
	std::cout<<"resized vectors results_en and results_ang to hold "<<((int)((cm_stop - cm_start)/cm_step))+1<<" elements"<<std::endl;
#endif

	for(int j=0;j<((int)(cm_stop - cm_start)/cm_step)+1;j++) {
		results_en.at(j) = -1;
		results_ang.at(j) = -1;
	}


	if(e_cm+q < ex)  {
		ex=e_cm+q;
		printf("# Readjusted excitation energy value to %f, since below threshold.\n",ex);
	}

	if(code == 0) {
#ifdef debug_orbits
		std::cout<<"#code 0: plot of lab energy as function of lab angles for ejectile"<<std::endl;
#endif
		for(cm=cm_start;cm<=cm_stop;cm +=cm_step) {
			lab2=lab2f(ex, cm, energy, particle);
			en2 = en2f(ex, cm, energy, particle);
			if(!(lab2>0) && !(lab2<0)) {
				lab2 = -1;
			}
#ifdef debug_orbits_results
			std::cout<<"lab2f("<<ex<<", "<<cm<<", "<<energy<<"): "<<lab2<<" "<<en2<<std::endl;
#endif
#ifdef debug_orbits
			std::cout<<ex<<"/"<<cm<<": results_en["<<(int)((cm - cm_start)/cm_step)<<"] = "<<flush;
#endif
			results_en.at((int)((cm - cm_start)/cm_step)) = en2;
			results_ang.at((int)((cm - cm_start)/cm_step)) = lab2/180.*TMath::Pi();
#ifdef debug_orbits
			std::cout<<en2<<", results_ang["<<(int)((cm - cm_start)/cm_step)<<"] = "<<lab2/180.*TMath::Pi()<<std::endl;
#endif
		}
	}
	else if(code == 1) {
#ifdef debug_orbits
		std::cout<<"#code 1: plot of lab energy as function of lab angles for recoil"<<std::endl;
#endif
		for(cm=cm_start;cm<=cm_stop;cm +=cm_step) {
			lab3=lab3f(ex, cm, energy, particle);
			en3 = en3f(ex, cm, energy, particle);
			if(!(lab3>0) && !(lab3<0)) {
				lab3 = -1;
			}
#ifdef debug_orbits_results
			std::cout<<"lab3f("<<ex<<", "<<cm<<", "<<energy<<"): "<<lab3<<" "<<en3<<std::endl;
#endif
#ifdef debug_orbits
			std::cout<<ex<<"/"<<cm<<": results_en["<<(int)((cm - cm_start)/cm_step)<<"] = "<<flush;
#endif
			results_en.at((int)((cm - cm_start)/cm_step)) = en3;
			results_ang.at((int)((cm - cm_start)/cm_step)) = lab3/180.*TMath::Pi();
#ifdef debug_orbits
			std::cout<<en3<<", results_ang["<<(int)((cm - cm_start)/cm_step)<<"] = "<<lab3/180.*TMath::Pi()<<std::endl;
#endif
		}
	}
	else if(code == 2) {
		std::cout<<"#code 2: plot of lab angles as function of cm angles"<<std::endl;
		for( cm = cm_start; cm <= cm_stop; cm += cm_step)  {
			lab2=lab2f(ex, cm, energy, particle); 
			printf("%6.2f %12.5e \n", cm, lab2);
		}
		for( cm = cm_start; cm <= cm_stop; cm += cm_step)  {
			lab3=lab3f(ex, cm, energy, particle); 
			printf("%6.2f %12.5e \n", cm, lab3);
		}
	}
	else if(code == 3) {
		// plot of lab energies as function of cm angles
		for( cm = cm_start; cm <= cm_stop; cm += cm_step)  {
			en2=en2f(ex, cm, energy, particle); 
			printf("%6.2f %12.5e \n", cm, en2);
		}
		for( cm = cm_start; cm <= cm_stop; cm += cm_step)  {
			en3=en3f(ex, cm, energy, particle); 
			printf("%6.2f %12.5e \n", cm, en3);
		}
	}

	return 0;
}


// included by S.Klupp
int Kinematic::orbits(Element* projectile, Element* target, Element* recoil, Element* ejectile, double energy, double cm, double ex, int code, double& results_en, double& results_ang) {
	double q,dq,e_cm,Ethr;
	double lab2,lab3,en2,en3;

	Element* particle[4] = {target, projectile, ejectile, recoil};

	Ethr=0;
	//load_table();  // load mass table

	CalculateQ(particle,&q,&dq);

	//
	e_cm = EnergyCm(energy,particle) - projectile->Mass() - target->Mass();

	// check for threshold
	if (e_cm < -q) {
		printf("# Energy below threshold\n");
		exit(1);
	}
	else if (q < -1.e-10) {
		// find threshold
		double hEtemp,lEtemp,e_cmt;
		hEtemp=energy;
		lEtemp=0;
		do {
			Ethr=(hEtemp+lEtemp)/2;
			e_cmt = EnergyCm(Ethr,particle) - projectile->Mass() - target->Mass();
			if (e_cmt < -q) {
				lEtemp=Ethr;
			}
			else {
				hEtemp=Ethr;
			}
		} while ( fabs(e_cmt + q)/(-q) > 1e-5);
		//printf("# Threshold = %f,  thr/A = %f \n",Ethr,Ethr/(projectile->Z() + projectile->N()));
		//std::cout << "q-value = " << q << std::endl;
	}

#ifdef debug_orbits
	printf("# Beam Energy = %f MeV\n",energy);
#endif

	// output

#ifdef debug_orbits
	std::cout<<"created new arrays results_en and results_ang with dimensions ["<<((int)((cm_stop - cm_start)/cm_step))+1<<"]"<<std::endl;
#endif

	results_en = -1;
	results_ang = -1;

	if(e_cm+q < ex) {
		ex=e_cm+q;
		printf("# Readjusted excitation energy value to %f, since below threshold.\n",ex);
	}

	if(code == 0) {
#ifdef debug_orbits
		std::cout<<"#code 0: plot of lab energy as function of lab angles for ejectile"<<std::endl;
#endif

		lab2=lab2f(ex, cm, energy, particle);
		en2 = en2f(ex, cm, energy, particle);
		if(!(lab2>0) && !(lab2<0)) {
			lab2 = -1;
		}

		results_en = en2;
		results_ang = lab2/180.*TMath::Pi();
	}
	else if(code == 1) {
#ifdef debug_orbits
		std::cout<<"#code 1: plot of lab energy as function of lab angles for recoil"<<std::endl;
#endif

		lab3=lab3f(ex, cm, energy, particle);
		en3 = en3f(ex, cm, energy, particle);
		if(!(lab3>0) && !(lab3<0)) {
			lab3 = -1;
		}

		results_en = en3;
		results_ang = lab3/180.*TMath::Pi();

	}
	else if(code == 2) {
		std::cout<<"#code 2: plot of lab angles as function of cm angles"<<std::endl;
		lab2=lab2f(ex, cm, energy, particle);
		printf("%6.2f %12.5e \n", cm, lab2);

		lab3=lab3f(ex, cm, energy, particle);
		printf("%6.2f %12.5e \n", cm, lab3);
	}
	else if(code == 3) {
		// plot of lab energies as function of cm angles
		en2=en2f(ex, cm, energy, particle);
		printf("%6.2f %12.5e \n", cm, en2);
		en3=en3f(ex, cm, energy, particle);
		printf("%6.2f %12.5e \n", cm, en3);
	}

	return 0;
}

//routine to get lowest cm-angle for given lab-angle (code == 0 => lab-angle of ejectile, code == 1 => lab-angle of recoil)
double Kinematic::CMAngle(double LabAngle, Element* projectile, Element* target, Element* recoil, Element* ejectile, int energies, double* energy, double ex, int code, double limit) {
	double cm1, cm2, lab1, lab2, lab3;
	int i=0;
	int counter;
	double result;
	double* cm;
	cm = new double[energies];

	Element* particle[4] = {target, projectile, ejectile, recoil};

#ifdef cmangle_debug
	std::cout<<target->A()<<target->Name()<<"("<<projectile->A()<<projectile->Name()<<","<<ejectile->A()<<ejectile->Name()<<")"<<recoil->A()<<recoil->Name()<<std::endl;
	std::cout<<energies<<std::endl;
	for(i=0;i<energies;i++) {
		std::cout<<energy[i]<<" ";
	}
	std::cout<<std::endl;
#endif
	//ejectile
	if(code == 0) {
		//std::cout << " ejectile " << std::endl;
		for(i=0;i<energies;i++) {
			//determine cm angle for ex
			counter = 0;
			cm1 = 0.01;
			cm2 = 179.99;
			lab1 = lab2f(ex, cm1, energy[i], particle)-LabAngle;
			lab2 = lab2f(ex, cm2, energy[i], particle)-LabAngle;
			if(lab1*lab2 >= 0) {
				lab3 = lab2f(ex, (cm1+cm2)/2., energy[i], particle)-LabAngle;
				if(lab1*lab3 < 0.) {
					cm2 = (cm1+cm2)/2.;
					lab2 = lab3;
				}
				else if(lab2*lab3 < 0.) {
					cm1 = (cm1+cm2)/2.;
					lab1 = lab3;
				}
				else {
					std::cerr<<"Error in Kinematic::CMAngle("<<LabAngle<<", "<<target->A()<<target->Name()<<", "<<projectile->A()<<projectile->Name()<<", "<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<", "<<energy[i]<<", "<<ex<<", "<<code<<", "<<limit<<"), LabAngle not bracketed by th_cm = "<<cm1<<" ("<<lab1+LabAngle<<"), th_cm = "<<cm2<<" ("<<lab2+LabAngle<<") or th_cm = "<<(cm1+cm2)/2.<<" ("<<lab3+LabAngle<<")"<<std::endl;
					return -1;
				}
			}
			while(fabs(lab1 - lab2) > limit) {
				lab3 = lab2f(ex, (cm1+cm2)/2., energy[i], particle)-LabAngle;
				if(lab1*lab3 < 0.) {
					cm2 = (cm1+cm2)/2.;
					lab2 = lab3;
				}
				else if(lab2*lab3 < 0.) {
					cm1 = (cm1+cm2)/2.;
					lab1 = lab3;
				}
				else {
					std::cerr<<"Error in Kinematic::CMAngle("<<LabAngle<<", "<<target->A()<<target->Name()<<", "<<projectile->A()<<projectile->Name()<<", "<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<", "<<energy[i]<<", "<<ex<<", "<<code<<", "<<limit<<"), LabAngle not bracketed by th_cm = "<<cm1<<" ("<<lab1+LabAngle<<"), th_cm = "<<cm2<<" ("<<lab2+LabAngle<<") or th_cm = "<<(cm1+cm2)/2.<<" ("<<lab3+LabAngle<<")"<<std::endl;
					return -1;
				}
				counter++;
				if(counter > 1000) {
					std::cerr<<"Maximum number of iterations reached"<<std::endl;
					break;
				}
#ifdef cmangle_debug
				std::cout<<lab1<<", "<<lab2<<", "<<lab3<<"; "<<cm1<<", "<<cm2<<", "<<(cm1+cm2)/2.<<std::endl;
#endif
			}

			cm[i] = (cm1+cm2)/2.;
		}
	}
	//recoil
	else if(code == 1) {
		//std::cout << " recoil " << std::endl;
		for(i=0;i<energies;i++) {
			//determine cm angle for ex
			counter = 0;
			cm1 = 0.01;
			cm2 = 179.99;
			//std::cout << "lab3f " << lab3f(ex, cm1, energy[i], particle)-LabAngle << " lab2f " << lab2f(ex, cm1, energy[i], particle)-LabAngle << std::endl;
			lab1 = lab3f(ex, cm1, energy[i], particle)-LabAngle;
			lab2 = lab3f(ex, cm2, energy[i], particle)-LabAngle;
			if(lab1*lab2 >= 0) {
				lab3 = lab2f(ex, (cm1+cm2)/2., energy[i], particle)-LabAngle;
				if(lab1*lab3 < 0.) {
					cm2 = (cm1+cm2)/2.;
					lab2 = lab3;
				}
				else if(lab2*lab3 < 0.) {
					cm1 = (cm1+cm2)/2.;
					lab1 = lab3;
				}
				else {
					std::cerr<<"Error in Kinematic::CMAngle("<<LabAngle<<", "<<target->A()<<target->Name()<<", "<<projectile->A()<<projectile->Name()<<", "<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<", "<<energy[i]<<", "<<ex<<", "<<code<<", "<<limit<<"), LabAngle not bracketed by th_cm = "<<cm1<<" ("<<lab1+LabAngle<<"), th_cm = "<<cm2<<" ("<<lab2+LabAngle<<") or th_cm = "<<(cm1+cm2)/2.<<" ("<<lab3+LabAngle<<")"<<std::endl;
					return -1;
				}
			}
			while(fabs(lab1 - lab2) > limit) {
				lab3 = lab3f(ex, (cm1+cm2)/2., energy[i], particle)-LabAngle;
				if(lab1*lab3 < 0.) {
					cm2 = (cm1+cm2)/2.;
					lab2 = lab3;
				}
				else if(lab2*lab3 < 0.) {
					cm1 = (cm1+cm2)/2.;
					lab1 = lab3;
				}
				else {
					std::cerr<<"Error in Kinematic::CMAngle("<<LabAngle<<", "<<target->A()<<target->Name()<<", "<<projectile->A()<<projectile->Name()<<", "<<ejectile->A()<<ejectile->Name()<<", "<<recoil->A()<<recoil->Name()<<", "<<energy[i]<<", "<<ex<<", "<<code<<", "<<limit<<"), LabAngle not bracketed by th_cm = "<<cm1<<" ("<<lab1+LabAngle<<"), th_cm = "<<cm2<<" ("<<lab2+LabAngle<<") or th_cm = "<<(cm1+cm2)/2.<<" ("<<lab3+LabAngle<<")"<<std::endl;
					return -1;
				}
				counter++;
				if(counter > 1000) {
					std::cerr<<"Maximum number of iterations reached"<<std::endl;
					break;
				}
#ifdef cmangle_debug
				std::cout<<lab1<<", "<<lab2<<", "<<lab3<<"; "<<cm1<<", "<<cm2<<", "<<(cm1+cm2)/2.<<std::endl;
#endif
			}

			cm[i] = (cm1+cm2)/2.;
		}
	}
	else {
		std::cerr<<"Error in Kinematic::CMAngle("<<LabAngle<<", "<<energy[i]<<", "<<ex<<", "<<code<<", "<<limit<<"), code can only be 0 or 1"<<std::endl;
		delete[] cm;

		return -1;
	}

	for(i=1;i<energies;i++) {
		if(cm[i] < cm[0]) {
			cm[0] = cm[i];
		}
	}

	result = cm[0];

	delete[] cm;

	return result;
}


//********************************************************************   
// 
//      PROGRAM  ENELOSS  (ION RANGES IN MATTER AND RELATED STUFF)       
// 
//      H. ERNST, ARGONNE, 1981  
//	with significant modification to the graphs and stopping power
// 	routines by K. LESKO Jan. 1984.
//
//********************************************************************   
// 

// rewritten in C++ by V. Bildstein 2005

Kinematic::Kinematic() {
}

Kinematic::Kinematic(Element* pro, Material* tar) {
	fProjectile = pro;
	fTarget = tar;
#ifdef debug_irma
	std::cout<<"Created instance of class Kinematic with projectile: "<<fProjectile->A()<<fProjectile->Name()<<std::endl
		<<"and target "<<fTarget->Name()<<std::endl;
	for(int i=0;i<fTarget->NumberOfElements();i++) {
		std::cout<<"    "<<i<<".: "<<fTarget->A(i)<<fTarget->Name(i)<<" ("<<fTarget->Fraction(i)<<")"<<std::endl;
	}
#endif
}

Kinematic::Kinematic(Element* pro, Material* tar, double target_thick) {
	fProjectile = pro;
	fTarget = tar;
	fTargetThickness = target_thick;
#ifdef debug_irma
	std::cout<<"Created instance of class Kinematic with projectile: "<<fProjectile->A()<<fProjectile->Name()<<std::endl
		<<"and target "<<fTarget->Name()<<" (thickness: "<<target_thick<<")"<<std::endl;
	for(int i=0;i<fTarget->NumberOfElements();i++) {
		std::cout<<"    "<<i<<".: "<<fTarget->A(i)<<fTarget->Name(i)<<" ("<<fTarget->Fraction(i)<<")"<<std::endl;
	}
#endif
}

Kinematic::~Kinematic() {
}

void Kinematic::SetThickness(double target_thick) {
	fTargetThickness = target_thick;
}

void Kinematic::Projectile(Element* pro) {
	fProjectile = pro;
}

void Kinematic::Target(Material* tar) {
	fTarget = tar;
}

int Kinematic::Run(double energy, double*& results) {
	//declaring neede variables
	double range;
	double dE_dx;
	double en_after_target;
	double en_loss;
	double dE_dx_after_target;
	double en_straggling;
	double ang_straggling;
	double chargestate;

#ifdef debug_irma
	std::cout<<"Kinematic::run: "<<energy<<", "<<fTargetThickness<<std::endl;
#endif

	//  if(!(fTargetThickness > 0) && !(fTargetThickness < 0))
	//    {
	//      std::cout<<"ERROR, fTargetThickness ("<<fTargetThickness<<") is not set!"<<std::endl;
	//
	//      return 1;
	//    }

	//calculate range
	range = CompoundRange(fProjectile, fTarget, energy, -7);

	//calculate dE/dx
#ifdef debug_irma
	std::cout<<"calling CompoundStoppingPower("<<fProjectile->Name()<<", "<<fTarget->Name()<<", "<<energy<<")"<<std::endl;
#endif

	dE_dx = CompoundStoppingPower(fProjectile, fTarget, energy);

	//calculate energy loss
	en_after_target = EnergyLoss(fProjectile, fTarget, fTargetThickness, energy, -4);

#ifdef debug_irma
	std::cout<<"en_after_target: "<<en_after_target<<std::endl;
#endif

	en_loss = energy - en_after_target;

#ifdef debug_irma
	std::cout<<"en_loss: "<<en_loss<<std::endl;
#endif

	//calculate dE/dx after target
	dE_dx_after_target = CompoundStoppingPower(fProjectile, fTarget, en_after_target);

	//calculate energy straggling
	en_straggling = EnergyStraggling(fProjectile, fTarget, dE_dx, dE_dx_after_target, en_loss);

	//calculate angular straggling
	ang_straggling = AngularStraggling(fProjectile, fTarget, fTargetThickness, energy);

	//calculate average charge state
	chargestate = ChargeState(fProjectile, energy);

	results = new double[9];

	results[0] = energy;
	results[1] = range;
	results[2] = en_loss;
	results[3] = en_straggling;
	results[4] = ang_straggling;
	results[5] = dE_dx;
	results[6] = dE_dx_after_target;
	results[7] = en_after_target;
	results[8] = chargestate;

	//std::cout<<"en_loss: "<<en_loss<<", results[2]: "<<results[2]<<", results: "<<results<<std::endl;

	return 0;
}

TSpline3* Kinematic::RangeVsEnergy(double maximum_energy, double step_size) {
	//declaring neede variables
	int NumberOfSteps = (int)(maximum_energy/step_size)+1;

	double* range = new double[NumberOfSteps];
	double* energy = new double[NumberOfSteps];
	double en;
	int i;

#ifdef debug_irma
	std::cout<<"Kinematic::RangeVsEnergy: "<<maximum_energy<<", "<<step_size<<std::endl;

	char tmpname[256];
	int write_file = snprintf(tmpname,256,"RangeVsEnergy_%d_%d.dat",fProjectile->N(),fProjectile->Z());
	std::ofstream tmp;
	if(write_file < 0 || write_file >= 256) {
		std::cerr<<"Error, string for filename too long, skipping writing to RangeVsEnergy_"<<fProjectile->N()<<"_"<<fProjectile->Z()<<".dat (write_file = "<<write_file<<")"<<std::endl;
		write_file = -1;
	}
	else
		tmp.open(tmpname);
#endif

	//calculate range
	for(i=0;i<NumberOfSteps;i++) {
		en = i*step_size;
		energy[i] = en;
		range[i] = CompoundRange(fProjectile, fTarget, en, -7);
#ifdef debug_irma
		//if(write_file > 0) ....
			tmp<<energy[i]<<" "<<range[i]<<std::endl;
			//std::cout<<energy[i]<<"<-- energy and range --> "<<range[i]<<" i: "<<i<<std::endl;
#endif
	}

#ifdef debug_irma
	if(write_file > 0)
		tmp.close();
#endif

	/* create TSpline3 object */
	//TGraph*   g      = new TGraph((int)(maximum_energy/step_size), energy, range);
	TSpline3* spline = new TSpline3("RangeVsEnergy", energy, range, NumberOfSteps);

	//delete g;
	delete[] range;
	delete[] energy;

	return spline;
}

TSpline3* Kinematic::EnergyVsRange(double maximum_energy, double step_size) {
	//declaring needed variables
	int NumberOfSteps = (int)(maximum_energy/step_size)+1;

	double* range = new double[NumberOfSteps];
	double* energy = new double[NumberOfSteps];
	double en;
	int i;

#ifdef debug_irma
	std::cout<<"Kinematic::EnergyVsRange: "<<maximum_energy<<", "<<step_size<<std::endl;

	char tmpname[256];
	int write_file = snprintf(tmpname,256,"EnergyVsRange_%d_%d.dat",fProjectile->N(),fProjectile->Z());
	std::ofstream tmp;
	if(write_file < 0 || write_file >= 256) {
		std::cerr<<"Error, string for filename too long, skipping writing to EnergyVsRange_"<<fProjectile->N()<<"_"<<fProjectile->Z()<<".dat (write_file = "<<write_file<<")"<<std::endl;
		write_file = -1;
	}
	else
		tmp.open(tmpname);
#endif

	//calculate range
	for(i=0;i<NumberOfSteps;i++) {
		en = i*step_size;
		energy[i] = en;
		range[i] = CompoundRange(fProjectile, fTarget, en, -7);
#ifdef debug_irma
		if(write_file > 0)
			tmp<<energy[i]<<" "<<range[i]<<std::endl;
#endif
	}

#ifdef debug_irma
	std::cout<<"test"<<std::endl;
	if(write_file > 0) {
		std::cout<<"closing tmp file"<<std::endl;
		tmp.close();
	}
	std::cout<<"creating graph"<<std::endl;
#endif

	/* create TSpline3 object */
	//TGraph*   g      = new TGraph((int)(maximum_energy/step_size), range, energy);
	TSpline3* spline = new TSpline3("EnergyVsRange", range, energy, NumberOfSteps);

#ifdef debug_irma
	std::cout<<"deleting stuff"<<std::endl;
#endif

	//possible without destroying TSpline ???
	//delete g;
	delete[] range;
	delete[] energy;

#ifdef debug_irma
	std::cout<<"return"<<std::endl;
#endif

	return spline;
}

TSpline3* Kinematic::EnergyVsThickness(double beam_energy, double step_size) {
#ifdef debug_VsX
	std::cout<<"\n\n\n\n Kinematic::EnergyVsThickness: "<<beam_energy<<", "<<step_size<<" target thickness: "<<fTargetThickness<<std::endl;
#endif

	//  if(!(fTargetThickness > 0) && !(fTargetThickness < 0))
	//    {
	//      std::cout<<"ERROR, fTargetThickness ("<<fTargetThickness<<") is not set!"<<std::endl;
	//
	//      return NULL;
	//    }

	//declaring needed variables
	int NumberOfSteps = (int)(fTargetThickness/step_size)+1;
	double* thickness;
	double* energy;
	if(NumberOfSteps > 1) {
		thickness = new double[NumberOfSteps];
		energy = new double[NumberOfSteps];
	}
	else {
		thickness = new double[2];
		energy = new double[2];
	}
	double thick;
	int i;
	//TGraph*   g = NULL;
	TSpline3* spline = NULL;

	//calculate range
	for(i=0;i<NumberOfSteps;i++) {
		thick = i*step_size;
		thickness[i] = thick;
		energy[i] = EnergyLoss(fProjectile, fTarget, thick, beam_energy, -4);
#ifdef debug_VsX
		std::cout<<"\n\n\n\n i: "<<i<<" thick: "<<thick<<" -> energy: "<<energy[i]<<" thikness/density: "<<thick/fTargetThickness<<std::endl;
#endif
	}

	if(NumberOfSteps == 0) {
		thickness[1] = thickness[0];
		energy[1] = energy[0];
		NumberOfSteps++;
	}

#ifdef debug_VsX
	std::cout<<"EnergyVsThickness: done with calculation of energy loss in "<<fTargetThickness<<" mg/cm2 and "<<(int)(fTargetThickness/step_size)<<" steps"<<std::endl;
#endif

	/* create TSpline3 object */
	//g      = new TGraph((int)(fTargetThickness/step_size), thickness, energy);
	spline = new TSpline3("EnergyVsThickness", thickness, energy, NumberOfSteps);

#ifdef debug_VsX
	std::cout<<"EnergyVsThickness: done with graph and spline"<<std::endl;
#endif

	//possible without destroying TSpline ???
	//delete g;
	delete[] thickness;
	delete[] energy;

#ifdef debug_VsX
	std::cout<<"EnergyVsThickness: done"<<std::endl;
#endif

	return spline;
}

TSpline3* Kinematic::RutherfordVsThickness(double beam_energy, double step_size, bool normalize) {
#ifdef debug_VsX
	std::cout<<"Kinematic::EnergyRange: "<<beam_energy<<", "<<step_size<<std::endl;
#endif

	//  if(!(fTargetThickness > 0) && !(fTargetThickness < 0))
	//    {
	//      std::cerr<<"ERROR, fTargetThickness ("<<fTargetThickness<<") is not set!"<<std::endl;
	//
	//      return NULL;
	//    }
	if(16*step_size > fTargetThickness && normalize) {
		std::cerr<<"Error, step_size is to big to normalize, setting step_size to "<<fTargetThickness/16.<<std::endl;
		step_size = fTargetThickness/16.;
	}

	//declaring needed variables
	int NumberOfSteps = (int)(fTargetThickness/step_size)+1;

#ifdef debug
	std::cout<<"allocating thickness and rutherford as double arrays of length "<<(int)(fTargetThickness/step_size)+1<<flush;
#endif

	double* thickness;
	double* rutherford;
	if(NumberOfSteps > 1) {
		thickness = new double[NumberOfSteps];
		rutherford = new double[NumberOfSteps];
	}
	else {
		thickness = new double[2];
		rutherford = new double[2];
	}

#ifdef debug
	std::cout<<" done"<<std::endl;
#endif
	double thick;
	double en;
	int i;

	//calculate range
	for(i=0;i<NumberOfSteps;i++) {
		thick = i*step_size;
		thickness[i] = thick;
		en = EnergyLoss(fProjectile, fTarget, thick, beam_energy, -4);
#ifdef debug_vsx
		std::cout<<"en ("<<i<<") = "<<en<<", pow(en,2) = "<<pow(en,2)<<", 1/pow(en,2) = "<<1/pow(en,2)<<std::endl;
#endif
		rutherford[i] = 1/pow(en,2);
	}

	if(NumberOfSteps == 0) {
		thickness[1] = thickness[0];
		rutherford[1] = rutherford[0];
		NumberOfSteps++;
	}

	//normalize if necessary
	if(normalize) {
		double integral[5];
		double h[5];
		for(i=0;i<5;i++) {
			integral[i] = 0.5*pow(rutherford[0],2);
			h[i] = pow(0.25,i);
#ifdef debug_vsx
			std::cout<<"integral["<<i<<"]: "<<integral[i]<<std::endl;
#endif
		}
		for(i=1;i<NumberOfSteps-1;i++) {
			if(i%16 == 0) {
				integral[0] += rutherford[i];
			}
			if(i%8 == 0) {
				integral[1] += rutherford[i];
			}
			if(i%4 == 0) {
				integral[2] += rutherford[i];
			}
			if(i%2 == 0) {
				integral[3] += rutherford[i];
			}
			integral[4] += rutherford[i];
		}
		for(i=0;i<5;i++) {
#ifdef debug_vsx
			std::cout<<"integral["<<i<<"]: "<<integral[i]<<std::endl;
#endif
			integral[i] += 0.5*rutherford[NumberOfSteps-1];
#ifdef debug_vsx
			std::cout<<"integral["<<i<<"]: "<<integral[i]<<std::endl;
#endif
			integral[i] *= step_size;
#ifdef debug_vsx
			std::cout<<"integral["<<i<<"]: "<<integral[i]<<std::endl;
#endif
		}

		int m,ns=4;
		double norm,ho,hp,w,den;
		double c[5],d[5];

		for(i=0;i<5;i++) {
			c[i]=integral[i];
			d[i]=integral[i];
		}
		norm=integral[ns--];
		for(m=1;m<5;m++)  { //For each column of the tableau,
			for(i=0;i<=4-m;i++)  { //we loop over the current c�s and d�s and update them. 
				ho=h[i];
				hp=h[i+m];
				w=c[i+1]-d[i];
				if((den=ho-hp) == 0.) {
					std::cerr<<"Error in RutherfordVsThickness, ho ("<<ho<<") == hp ("<<hp<<")"<<std::endl;
					exit(1);
				}
				den=w/den;
				d[i]=hp*den; 
				c[i]=ho*den;
			}
			norm += d[ns--]; // d[ns--] is also dnorm
		}

#ifdef debug_vsx
		std::cout<<"RutherfordVsX: norm = "<<norm<<std::endl;
		for(i=0;i<5;i++) {
			std::cout<<"integral["<<i<<"] = "<<integral[i]<<std::endl;
		}
#endif

		for(i=0;i<NumberOfSteps;i++) {
#ifdef debug_vsx
			std::cout<<"change rutherford["<<i<<"] from "<<rutherford[i]<<" to "<<rutherford[i]/norm<<std::endl;
#endif
			rutherford[i] /= norm;
		}
	} //if(normalize)

	/* create TSpline3 object */
	//TGraph*   g      = new TGraph(NumberOfSteps, thickness, rutherford);
	TSpline3* spline = new TSpline3("RutherfordVsThickness", thickness, rutherford, NumberOfSteps);

	//possible without destroying TSpline ???
	//delete g;
	delete[] rutherford;
	delete[] thickness;

	return spline;
}

TSpline3* Kinematic::EnergyLossVsEnergy(double max_energy, double step_size, double factor) {
#ifdef debug_VsE
	std::cout<<"Kinematic::EnergyLossVsEnergy: "<<max_energy<<", "<<step_size<<std::endl;
	std::cout<<"projectile = "<<fProjectile->Name()<<", A = "<<fProjectile->A()<<", N = "<<fProjectile->N()<<", Z = "<<fProjectile->Z()<<std::endl
		<<"target = "<<fTarget->Name()<<", thickness = "<<fTargetThickness<<std::endl;
#endif

	TSpline3* spline = NULL;

	//  if(!(fTargetThickness > 0) && !(fTargetThickness < 0))
	//    {
	//      std::cout<<"ERROR, fTargetThickness ("<<fTargetThickness<<") is not set!"<<std::endl;
	//
	//      return spline;
	//    }

	//declaring needed variables
	int NumberOfSteps = (int)(max_energy/step_size)+1;

#ifdef debug_VsE
	std::cout<<"arrays of size "<<NumberOfSteps<<std::endl;
#endif

	double* energyloss;
	double* energy;
	if(NumberOfSteps > 1) {
		energyloss = new double[NumberOfSteps];
		energy = new double[NumberOfSteps];
	}
	else {
		energyloss = new double[2];
		energy = new double[2];
	}

	int i;

	//calculate energy loss
	for(i = 0;i < NumberOfSteps; i++) {
		energy[i] = i*step_size*factor;
		energyloss[i] = EnergyLoss(energy[i]/factor)*factor;
	}

	if(NumberOfSteps == 0) {
		energyloss[1] = energyloss[0];
		energy[1] = energy[0];
		NumberOfSteps++;
	}

#ifdef debug_VsE
	std::cout<<"EnergyLossVsEnergy: done with calculation of energy loss in "<<fTargetThickness<<" mg/cm2 and "<<NumberOfSteps<<" steps"<<std::endl;
#endif

	/* create TSpline3 object */
	spline = new TSpline3(Form("EnergyLossVsEnergy_%s_%s",fProjectile->Name(),fTarget->Name()), energy, energyloss, NumberOfSteps);

#ifdef debug_VsE
	std::cout<<"EnergyLossVsEnergy: done with spline"<<std::endl;
#endif

	//possible without destroying TSpline ???
	delete[] energy;
	delete[] energyloss;

#ifdef debug_VsE
	std::cout<<"EnergyLossVsEnergy: done"<<std::endl;
#endif

	return spline;
}

void Kinematic::EnergyLossVsEnergy(double max_energy, std::vector<double>& energy, std::vector<double>& energyloss, double step_size, double factor) {
#ifdef debug_VsE
	std::cout<<"Kinematic::EnergyLossVsEnergy: "<<max_energy<<", "<<step_size<<std::endl;
	std::cout<<"projectile = "<<fProjectile->Name()<<", A = "<<fProjectile->A()<<", N = "<<fProjectile->N()<<", Z = "<<fProjectile->Z()<<std::endl
		<<"target = "<<fTarget->Name()<<", thickness = "<<fTargetThickness<<std::endl;
#endif

	//  if(!(fTargetThickness > 0) && !(fTargetThickness < 0))
	//    {
	//      std::cout<<"ERROR, fTargetThickness ("<<fTargetThickness<<") is not set!"<<std::endl;
	//
	//      return;
	//    }

	//declaring needed variables
	int NumberOfSteps = (int)(max_energy/step_size)+1;

#ifdef debug_VsE
	std::cout<<"arrays of size "<<NumberOfSteps<<std::endl;
#endif

	energy.resize(NumberOfSteps);
	energyloss.resize(NumberOfSteps);

	unsigned int i;

	//calculate energy loss
	for(i = 0; i < energy.size(); i++) {
		energy.at(i) = i*step_size*factor;
		energyloss.at(i) = EnergyLoss(energy.at(i)/factor)*factor;
	}

#ifdef debug_VsE
	std::cout<<"EnergyLossVsEnergy: done with calculation of energy loss in "<<fTargetThickness<<" mg/cm2 and "<<(int)(max_energy/step_size)<<" steps"<<std::endl;
#endif

	return;
}

double Kinematic::Range(double energy) {
	return CompoundRange(fProjectile, fTarget, energy, -7);
}

double Kinematic::EnergyLoss(double energy) {
	return (energy - EnergyLoss(fProjectile, fTarget, fTargetThickness, energy, -4));
}

double Kinematic::EnergyAfterTarget(double energy) {
	return EnergyLoss(fProjectile, fTarget, fTargetThickness, energy, -4);
}

double Kinematic::Energy(double range, double max_energy, double tolerance) {
	//std::cout<<std::endl<<"Energy("<<range<<", "<<max_energy<<", "<<tolerance<<")"<<std::endl;
	double tmp_range, tmp_energy;

	tmp_energy = max_energy;
	tmp_range = CompoundRange(fProjectile, fTarget, tmp_energy, -7);

	if(tmp_range < range) {
		return -1.;
	}

	while(fabs(tmp_range-range) > tolerance) {
		tmp_energy *= range/tmp_range;
		tmp_range = CompoundRange(fProjectile, fTarget, tmp_energy, -7);
	}

	//std::cout<<"done: "<<tmp_range<<" - "<<range<<" = "<<fabs(tmp_range-range)<<" <= "<<tolerance<<std::endl;

	return tmp_energy;
}

//double Kinematic::ReverseEnergy(double DeltaE, double E, double theta, double*& energies)
//{
//  return 0.;
//}

double Kinematic::CompoundRange(Element* projectile, Material* target, double energy, int integration_limit) {
#ifdef debug_irma
	std::cout<<"starting CompoundRange("<<std::flush<<projectile->A()<<projectile->Name()<<", "<<std::flush<<target->Name()<<", "<<energy<<", "<<integration_limit<<")"<<std::endl;
#endif
	int z_p = projectile->Z();
	double a_p = projectile->A();
	double en_1;
	double en_2;
	int n;
	double range;
	double h;

	en_1 = pow(10.,integration_limit);

#ifdef debug_irma
	std::cout<<"en_1: "<<en_1<<std::endl;
#endif

	if(energy < en_1 || z_p <= 0 || a_p <= 0) {
		return 0;
	}

	n = -integration_limit - 1 + (int) floor(log(energy)*0.4342944819 + 0.5);
	en_2 = energy/pow(10,(double) n);

#ifdef debug_irma
	std::cout<<"en_2: "<<en_2<<std::endl;
#endif

	//integrate using Weddle formula
	range = 1e-30;
	n++;

#ifdef debug_irma
	std::cout<<"starting loop from 0 to n-1 = "<<n-1<<std::endl;
#endif
	for(int i=0;i<n;i++) {
		h = (en_2 - en_1)/6.;
		range += h*(0.3/CompoundStoppingPower(projectile, target, en_1)       + 1.5/CompoundStoppingPower(projectile, target, en_1 + h)
				+ 0.3/CompoundStoppingPower(projectile, target, en_1 + 2*h) + 1.8/CompoundStoppingPower(projectile, target, en_1 + 3*h)
				+ 0.3/CompoundStoppingPower(projectile, target, en_2 - 2*h) + 1.5/CompoundStoppingPower(projectile, target, en_2 - h)
				+ 0.3/CompoundStoppingPower(projectile, target, en_2));
		en_1 = en_2;
		en_2 *= 10;
#ifdef debug_irma
		std::cout<<"Done with iteration "<<i<<", range = "<<range<<std::endl;
#endif
	}
#ifdef debug_irma
	std::cout<<"Done with loop, returning CompoundRange = "<<range<<std::endl;
#endif

	return range;
}

double Kinematic::CompoundStoppingPower(Element* projectile, Material* target, double energy) {
#ifdef debug_irma
	std::cout<<"starting CompoundStoppingPower("<<std::flush<<projectile->A()<<projectile->Name()<<", "<<std::flush<<target->Name()<<", "<<energy<<")"<<std::endl;
#endif
	double res;

	res = 0.;

	//std::cout<<"cstop: NumberOfElements: "<<target->NumberOfElements()<<std::endl;
	for(size_t i = 0; i < target->NumberOfElements(); ++i) {
		res += target->Fraction(i)*StoppingPower(projectile, target->GetElement(i), energy, target->Gaseous());
#ifdef debug_results
		std::cout<<"CSTOP: "<<target->NumberOfElements()<<" "<<res<<" "<<target->Fraction(i)<<std::endl;
#endif
	}

	if(res < 1.e-30) {
		res = 1.e-30;
	}

	return res;
}

double Kinematic::StoppingPower(Element* projectile, Element* target, double energy, bool gaseous) {
#ifdef debug_irma
	std::cout<<"starting StoppingPower("<<projectile->A()<<projectile->Name()<<", "<<target->A()<<target->Name()<<", "<<energy<<", "<<gaseous<<")"<<std::endl;
#endif
	int z_p = projectile->Z();
	double a_p = projectile->A();
	int z_t = target->Z();
	double a_t = target->A();

	if(z_p <= 0 || a_p <= 0 || z_t <= 0 || a_t <= 0 || energy <= 1.e-30) {
		return 0;
	}

	if(z_t > 92) {
		z_t = 92;
	}

	double vsq;
	double vovc;
	double factor;  //to convert units from eV/10^15 Atoms/cm^2 to MeV/mg/cm^2
	double a;
	double eps;
	double dedx_n;
	double dedx_e;
	double e_h;
	//double bsq;
	double z_eff;
	double hexp;
	double beta;
	double e_he;

	vsq = energy/a_p;
	vovc = 0.04634*sqrt(vsq);
	if(vovc > 1) {
		vovc = 0.999999;
	}
	factor = 0.60225/a_t;

	//nuclear stopping powers
	a = a_p + a_t;
	z_eff = pow(z_p,0.23) + pow(z_t,0.23);
	eps = 32520.*energy*a_t/(z_eff*z_p*z_t*a);
	if(eps > 30.) {
		dedx_n = 0.5*log(eps)/eps;
	} else {
		dedx_n = log(1. + 1.1383*eps)/(2.*(eps + 0.01321*pow(eps,0.21226) + 0.19593*sqrt(eps)));
#ifdef debug_results
		std::cout<<"STOP1: "<<std::setprecision(10)<<log(1. + 1.1383*eps)<<" "<<(1. + 1.1383*eps)<<" "<<eps<<" "<<0.01321*pow(eps,0.21226)<<" "<<0.19593*sqrt(eps)<<std::endl;
#endif
	}
#ifdef debug_results
	std::cout<<"STOP2: "<<z_p<<" "<<a_p<<" "<<z_t<<" "<<a_t<<std::setprecision(10)<<" "<<energy<<" "<<vsq<<" "<<vovc<<" "<<dedx_n<<" "<<eps<<std::endl;
#endif
	dedx_n *= 8.462*z_p*z_t*a_p/(a*z_eff)*factor;
#ifdef debug_results
	std::cout<<"STOP3: "<<std::setprecision(10)<<dedx_n<<" "<<8.462*z_p*z_t*a_p/(a*z_eff)*factor<<std::endl;
#endif

	//master electronic stopping powers for hydrogen
	if(z_p != 2) {
		e_h = 1000.*vsq;

		// new version
		dedx_e = factor/(1./(HydrogenA(2,z_t)*pow(e_h,0.45)) + e_h/(HydrogenA(3,z_t)*log(1. + HydrogenA(4,z_t)/e_h + HydrogenA(5,z_t)*e_h)));

		//calculate Ziegler's effective Z's and electronic stopping powers (nothing to be done for hydrogen and special case Lithium)
		if(z_p == 1) {
			return (dedx_n + dedx_e);
		} else if(z_p == 3) {
			z_eff = 3.;
			hexp = 0.7138 + 0.002797*e_h + 1.348e-6*e_h*e_h;
			if(hexp < 60) {
				z_eff -= 3./exp(hexp);
			}
			dedx_e *= z_eff*z_eff;
#ifdef debug_results
			std::cout<<"STOP5.2: "<<dedx_e<<" "<<dedx_n<<std::endl;
#endif
			return (dedx_n + dedx_e);
		} else {
			beta = 0.1772*sqrt(e_h)/pow(z_p,2./3.);
			z_eff = z_p*(1. - (1.034 - 0.1777/exp(0.08114*z_p))/exp(beta + 0.0378*sin(1.5708*beta)));
			if(z_eff < 1.) {
				z_eff = 1.;
			}
			dedx_e *= z_eff*z_eff;
#ifdef debug_results
			std::cout<<"STOP5.1: "<<dedx_e<<" "<<dedx_n<<std::endl;
#endif
			return (dedx_n + dedx_e);
		}
	} else {
		e_he = 4.*vsq;
		if(e_he >= 10.) { //solid and gaseous matter , high energies
			dedx_e = factor*exp(HeliumA(6,z_t) + (HeliumA(7,z_t) + (HeliumA(8,z_t) + HeliumA(9,z_t)*log(1./e_he))*log(1./e_he))*log(1./e_he));
#ifdef debug_results
			std::cout<<"STOP5.3: "<<dedx_e<<" "<<dedx_n<<std::endl;
#endif
			return (dedx_n + dedx_e);
		} else { //low energies
			if(!gaseous) { //solid matter
				dedx_e = factor/(1./(pow(1000.*e_he,HeliumA(2,z_t))*HeliumA(1,z_t))+e_he/(HeliumA(3,z_t)*log(1.+HeliumA(4,z_t)/e_he+HeliumA(5,z_t)*e_he)));
#ifdef debug_results
				std::cout<<"STOP5.4: "<<dedx_e<<" "<<dedx_n<<std::endl;
#endif
				return (dedx_n + dedx_e);
			} else { //gaseous matter
				dedx_e = factor/(1./(pow(1000.*e_he,HeliumB(2,z_t))*HeliumB(1,z_t))+e_he/(HeliumB(3,z_t)*log(1.+HeliumB(4,z_t)/e_he+HeliumB(5,z_t)*e_he)));
#ifdef debug_results
				std::cout<<"STOP5.5: "<<dedx_e<<" "<<dedx_n<<std::endl;
#endif
				return (dedx_n + dedx_e);
			}
		}
	}
}

double Kinematic::EnergyLoss(Element* projectile, Material* target, double target_thickness, double energy, int integration_limit) {
#ifdef debug_irma
	std::cout<<"started EnergyLoss("<<projectile->A()<<projectile->Name()<<", "<<target->Name()<<", "<<target_thickness<<", "<<energy<<", "<<integration_limit<<")"<<std::endl;
#endif
	double en;
	double range;
	double dedx;
	double dx;
	double x;

	if(target_thickness < 0) {
		std::cout<<"Warning, target thickness < 0: EnergyLoss("<<projectile->A()<<projectile->Name()<<", "<<target->Name()<<", "<<target_thickness<<", "<<energy<<", "<<integration_limit<<"), returning input energy"<<std::endl;
		return energy;
	}

	if(target_thickness == 0) {
		//no target => no energy loss
		return energy;
	}

	en = energy;

	x = 0;

	if(energy <= 0) {
		range = 0.;
		dedx = 100000.;
	}
	else {
		range = CompoundRange(projectile, target, energy, integration_limit);
		dedx = CompoundStoppingPower(projectile, target, energy);
	}

	dx = 0.1/dedx;

	if((range - target_thickness) <= 0) {
		return 0.;
	}

	while(1) {
		if(target_thickness - x - dx < 0) {
#ifdef debug_irma
			std::cout<<"target_thickness: "<<target_thickness<<", x: "<<x<<", dx: "<<dx<<std::endl;
#endif
			return (en - dedx*(target_thickness - x));
		}
		en -= dedx*dx;
#ifdef debug_irma
		std::cout<<"target_thickness: "<<target_thickness<<", x: "<<x<<", dx: "<<dx<<", en: "<<en<<std::endl;
#endif
		x += dx;
		dedx = CompoundStoppingPower(projectile, target, en);
	}
}

double Kinematic::EnergyStraggling(Element* projectile, Material* target, double dE_dx, double dE_dx_after_target, double energy_loss) {
	//see H. Schmidt-Boeking in Lecture Notes in Physics 83 (1978)

	int z_p = projectile->Z();
	double a_p = projectile->A();

	int n_el = target->NumberOfElements();

	if(z_p < 0 || a_p <= 0 || n_el < 1 || n_el > 10 || dE_dx < 1.e-10) {
		return 0;
	}

	double z_t = 0;

	for(int i=0;i<n_el;i++) {
		z_t += target->Z(i)*target->Fraction(i);
	}

	double z_eff;

	z_eff = pow((double) z_p, 1./3.) + pow(z_t, 1./3.);

#ifdef e_strag_debug
	std::cout<<"ESTRAG: "<<sqrt(z_p*z_t/z_eff)<<" "<<pow(energy_loss,0.53)<<" "<<dE_dx_after_target/dE_dx<<std::endl;
#endif

	return (0.014*sqrt(z_p*z_t/z_eff)*pow(energy_loss,0.53)*dE_dx_after_target/dE_dx);
}

double Kinematic::AngularStraggling(Element* projectile, Material* target, double target_thickness, double energy) {
	//see L.Meyer, Phys. Stat. Sol. 44, 253(1971)

	double z_p = projectile->Z();
	double a_p = projectile->A();

	int n_el = target->NumberOfElements();

	if(z_p < 0 || a_p < 0 || n_el < 1 || n_el > 10 || energy < 1e-30) {
		return 0;
	}

	double ang = 0;
	double z_eff;
	double tau;
	double exponent;

	for(int i=0;i<n_el;i++) {
		z_eff = sqrt(pow(z_p, 2./3.)+pow(target->Z(i), 2./3.));
		tau = 41470.0*target_thickness/(pow(z_eff,2)*target->A(i));
		exponent = pow(log(1.03+tau),-0.115) - 0.115;
		ang += target->Fraction(i)*target->Z(i)*z_eff*pow(tau,exponent); // should be a += ??
	}
	if(energy < 1.e-30) {
		return (0.0022*z_p*ang/(energy+1.e-30));
	}
	else {
		return (0.0022*z_p*ang/energy);
	}
}

double Kinematic::ChargeState(Element* projectile, double energy) {
	//see H.-D. Betz, Rev. Mod. Phys 44, 465(1972)

	double z = projectile->Z();
	double a = projectile->A();
	if(energy < 0.0001 || a <= 0 || z <= 0) {
		return 0;
	}

	return z*pow(1.0+0.105*pow(z,0.75)*pow(a/energy,0.83),-0.6);
}

double Kinematic::ShellCorrection(int z) {
	//      SHELL CORRECTIONS FOR BETHE FORMULA
	//      FOR HYDROGEN AS PROJECTILE.        
	switch(z) {
		case 1:
			return 4.;
			break;
		case 2:
			return -3.;
			break;
		case 3:
			return -9.;
			break;
		case 4:
			return -8.;
			break;
		case 5:
			return -7.;
			break;
		case 6:
			return -9.;
			break;
		case 7:
			return -6.;
			break;
		case 8:
			return -9.;
			break;
		case 9:
			return -10.;     
			break;
		case 10:
			return -11.;    
			break;
		default:
			if(11 <= z && z <= 50) {
				return (-21.258 + 2.1282*z - 0.099659*z*z + 0.0011623*z*z*z + 0.000024822*z*z*z*z - 0.0000004405*z*z*z*z*z);
			}
			else if(51 <= z && z <= 92) {
				return (1089.1  - 85.686*z + 2.6457*z*z   - 0.040082*z*z*z  + 0.0002976*z*z*z*z   - 0.00000086478*z*z*z*z*z);
			}
			break;
	}
	return -1.;
}

double Kinematic::HydrogenA(int index, int z) {
	double res[92][7] = {{1.262, 1.440,  242.6,12000.0,0.115900,0.000510,54360.0},
								{1.229, 1.397,  484.5, 5873.0,0.052250,0.001020,24510.0},
								{1.411, 1.600,  725.6, 3013.0,0.045780,0.001530,21470.0},
								{2.248, 2.590,  966.0,  153.8,0.034750,0.002039,16300.0},
								{2.474, 2.815, 1206.0, 1060.0,0.028550,0.002549,13450.0},
								{2.631, 2.989, 1445.0,  957.2,0.028190,0.003059,13220.0},
								{2.954, 3.350, 1683.0, 1900.0,0.025130,0.003569,11790.0},
								{2.652, 3.000, 1920.0, 2000.0,0.022300,0.004079,10460.0},
								{2.085, 2.352, 2157.0, 2634.0,0.018160,0.004589, 8517.0},
								{1.951, 2.199, 2393.0, 2699.0,0.015680,0.005099, 7353.0},
								{2.542, 2.869, 2628.0, 1854.0,0.014720,0.005609, 6905.0},
								{3.792, 4.293, 2862.0, 1009.0,0.013970,0.006118, 6551.0},
								{4.154, 4.739, 2766.0,  164.5,0.020230,0.006628, 6309.0},
								{4.150, 4.700, 3329.0,  550.0,0.013210,0.007138, 6194.0},
								{3.232, 3.647, 3561.0, 1560.0,0.012670,0.007648, 5942.0},
								{3.447, 3.891, 3792.0, 1219.0,0.012110,0.008158, 5678.0},
								{5.047, 5.714, 4023.0,  878.6,0.011780,0.008668, 5524.0},
								{5.731, 6.500, 4253.0,  530.0,0.011230,0.009178, 5268.0},
								{5.151, 5.833, 4482.0,  545.7,0.011290,0.009687, 5295.0},
								{5.521, 6.252, 4710.0,  553.3,0.011120,0.010200, 5214.0},
								{5.201, 5.884, 4938.0,  560.9,0.009995,0.010710, 4688.0},
								{4.862, 5.496, 5165.0,  568.5,0.009474,0.011220, 4443.0},
								{4.480, 5.055, 5391.0,  952.3,0.009117,0.011730, 4276.0},
								{3.983, 4.489, 5616.0, 1336.0,0.008413,0.012240, 3946.0},
								{3.469, 3.907, 5725.0, 1461.0,0.008829,0.012750, 3785.0},
								{3.519, 3.963, 6065.0, 1243.0,0.007782,0.013260, 3650.0},
								{3.140, 3.535, 6288.0, 1372.0,0.007361,0.013770, 3453.0},
								{3.553, 4.004, 6205.0,  555.1,0.008763,0.014280, 3297.0},
								{3.696, 4.175, 4673.0,  387.8,0.021880,0.014790, 3174.0},
								{4.210, 4.750, 6953.0,  295.2,0.006809,0.015300, 3194.0},
								{5.041, 5.697, 7173.0,  202.6,0.006725,0.015810, 3154.0},
								{5.554, 6.300, 6496.0,  110.0,0.009689,0.016320, 3097.0},
								{5.323, 6.012, 7611.0,  292.5,0.006447,0.016830, 3024.0},
								{5.874, 6.656, 7395.0,  117.5,0.007684,0.017340, 3006.0},
								{5.611, 6.335, 8046.0,  365.2,0.006244,0.017850, 2928.0},
								{6.411, 7.250, 8262.0,  220.0,0.006087,0.018360, 2855.0},
								{5.694, 6.429, 8478.0,  292.9,0.006087,0.018860, 2855.0},
								{6.339, 7.159, 8693.0,  330.3,0.006003,0.019370, 2815.0},
								{6.407, 7.234, 8907.0,  367.8,0.005889,0.019880, 2762.0},
								{6.734, 7.603, 9120.0,  405.2,0.005765,0.020390, 2700.0},
								{6.902, 7.791, 9333.0,  442.7,0.005587,0.020900, 2621.0},
								{6.425, 7.248, 9545.0,  480.2,0.005367,0.021410, 2517.0},
								{6.799, 7.671, 9756.0,  517.6,0.005315,0.021920, 2493.0},
								{6.108, 6.887, 9966.0,  555.1,0.005151,0.022430, 2416.0},
								{5.924, 6.677,10180.0,  592.5,0.004919,0.022940, 2307.0},
								{5.238, 5.900,10380.0,  630.0,0.004758,0.023450, 2231.0},
								{5.623, 6.354, 7160.0,  337.6,0.013940,0.023960, 2193.0},
								{5.814, 6.554,10800.0,  355.5,0.004626,0.024470, 2170.0},
								{6.230, 7.024,11010.0,  370.9,0.004540,0.024980, 2129.0},
								{6.410, 7.227,11210.0,  386.4,0.004474,0.025490, 2099.0},
								{7.500, 8.480, 8608.0,  348.0,0.009074,0.026000, 2069.0},
								{6.979, 7.871,11620.0,  392.4,0.004402,0.026510, 2065.0},
								{7.725, 8.716,11830.0,  394.8,0.004376,0.027020, 2052.0},
								{8.231, 9.289,12030.0,  397.3,0.004384,0.027530, 2056.0},
								{7.287, 8.218,12230.0,  399.7,0.004447,0.028040, 2086.0},
								{7.899, 8.911,12430.0,  402.1,0.004511,0.028550, 2116.0},
								{8.041, 9.071,12630.0,  404.5,0.004540,0.029060, 2129.0},
								{7.489, 8.444,12830.0,  406.9,0.004420,0.029570, 2073.0},
								{7.291, 8.219,13030.0,  409.3,0.004298,0.030080, 2016.0},
								{7.098, 8.000,13230.0,  411.8,0.004182,0.030590, 1962.0},
								{6.910, 7.786,13430.0,  414.2,0.004058,0.031100, 1903.0},
								{6.728, 7.580,13620.0,  416.6,0.003976,0.031610, 1865.0},
								{6.551, 7.380,13820.0,  419.0,0.003877,0.032120, 1819.0},
								{6.739, 7.592,14020.0,  421.4,0.003863,0.032630, 1812.0},
								{6.212, 7.996,14210.0,  423.9,0.003725,0.033140, 1747.0},
								{5.517, 6.210,14440.0,  426.3,0.003632,0.033650, 1703.0},
								{5.219, 5.874,14600.0,  428.7,0.003498,0.034160, 1640.0},
								{5.071, 5.706,14790.0,  433.0,0.003405,0.034670, 1597.0},
								{4.926, 5.542,14980.0,  433.5,0.003342,0.035180, 1567.0},
								{4.787, 5.386,15170.0,  435.9,0.003292,0.035690, 1544.0},
								{4.893, 5.505,15360.0,  438.4,0.003243,0.036200, 1521.0},
								{5.028, 5.657,15550.0,  440.8,0.003195,0.036710, 1499.0},
								{4.738, 5.329,15740.0,  443.2,0.003186,0.037220, 1494.0},
								{4.574, 5.144,15930.0,  442.4,0.003144,0.037730, 1475.0},
								{5.200, 5.851,16120.0,  441.6,0.003122,0.038240, 1464.0},
								{5.070, 5.704,16300.0,  440.9,0.003082,0.038750, 1446.0},
								{4.945, 5.563,16490.0,  440.1,0.002965,0.039260, 1390.0},
								{4.476, 5.034,16670.0,  439.3,0.002871,0.039770, 1347.0},
								{4.856, 5.460,18320.0,  438.5,0.002542,0.040280, 1354.0},
								{4.308, 4.843,17040.0,  487.8,0.002882,0.040790, 1352.0},
								{4.723, 5.311,17220.0,  537.0,0.002913,0.041300, 1366.0},
								{5.319, 5.982,17400.0,  586.3,0.002871,0.041810, 1347.0},
								{5.956, 6.700,17800.0,  677.0,0.002660,0.042320, 1336.0},
								{6.158, 6.928,17770.0,  586.3,0.002813,0.042830, 1319.0},
								{6.204, 6.979,17950.0,  586.3,0.002776,0.043340, 1302.0},
								{6.181, 6.954,18120.0,  586.3,0.002748,0.043850, 1289.0},
								{6.949, 7.820,18300.0,  586.3,0.002737,0.044360, 1284.0},
								{7.506, 8.448,18480.0,  586.3,0.002727,0.044870, 1279.0},
								{7.649, 8.609,18660.0,  586.3,0.002697,0.045380, 1265.0},
								{7.710, 8.679,18830.0,  586.3,0.002641,0.045890, 1239.0},
								{7.407, 8.336,19010.0,  586.3,0.002603,0.046400, 1221.0},
								{7.290, 8.204,19180.0,  586.3,0.002573,0.046910, 1207.0}};

	return res[z-1][index-1];
}

double Kinematic::HeliumA(int index, int z) {
	double res [92][9] = {{ 0.9661,  0.4126,  6.9200,  8.8310,  2.5820, 2.371000, 0.546200,-0.079320,-0.006853},
								 { 2.0270,  0.2931, 26.3400,  6.6600,  0.3409, 2.809000, 0.484700,-0.087560,-0.007281}, 
								 { 1.4200,  0.4900, 12.2500, 32.0000,  9.1610, 3.095000, 0.443400,-0.092590,-0.007459}, 
								 { 2.2060,  0.5100, 15.3200,  0.2500,  8.9950, 3.280000, 0.418800,-0.095640,-0.007604}, 
								 { 3.6910,  0.4128, 18.4800, 50.7200,  9.0000, 3.426000, 0.400000,-0.097960,-0.007715},
								 { 4.2320,  0.3877, 22.9900, 35.0000,  7.9930, 3.588000, 0.392100,-0.099350,-0.007804},
								 { 2.5100,  0.4752, 38.2600, 13.0200,  1.8920, 3.759000, 0.409400,-0.096460,-0.007661},
								 { 1.7660,  0.5261, 37.1100, 15.2400,  2.8040, 3.782000, 0.373400,-0.101100,-0.007874},
								 { 1.5330,  0.5310, 40.4400, 18.4100,  2.7180, 3.816000, 0.350400,-0.104600,-0.008074},
								 { 1.1830,  0.5500, 39.8300, 17.4900,  4.0010, 3.863000, 0.334200,-0.107200,-0.008231},
								 { 9.8940,  0.3081, 23.6500,  0.3840, 92.9300, 3.898000, 0.319100,-0.108600,-0.008271},
								 { 4.3000,  0.4700, 34.3000,  3.3000, 12.7400, 3.961000, 0.314000,-0.109100,-0.008297},
								 { 2.5000,  0.6250, 45.7000,  0.1000,  4.3590, 4.024000, 0.311300,-0.109300,-0.008306},
								 { 2.1000,  0.6500, 49.3400,  1.7880,  4.1330, 4.077000, 0.307400,-0.108900,-0.008219},
								 { 1.7290,  0.6562, 53.4100,  2.4050,  3.8450, 4.124000, 0.302300,-0.109400,-0.008240},
								 { 1.4020,  0.6791, 58.9800,  3.5280,  3.2110, 4.164000, 0.296400,-0.110100,-0.008267},
								 { 1.1170,  0.7044, 69.6900,  3.7050,  2.1560, 4.210000, 0.293600,-0.110300,-0.008270},
								 { 0.9172,  0.7240, 79.4400,  3.6480,  1.6460, 4.261000, 0.299400,-0.108500,-0.008145},
								 { 8.5540,  0.3817, 83.6100, 11.8400,  1.8750, 4.300000, 0.290300,-0.110300,-0.008259},
								 { 6.2970,  0.4622, 65.3900, 10.1400,  5.0360, 4.344000, 0.289700,-0.110200,-0.008245},
								 { 5.3070,  0.4918, 61.7400, 12.4000,  6.6650, 4.327000, 0.270700,-0.112700,-0.008370},
								 { 4.7100,  0.5087, 65.2800,  8.8060,  5.9480, 4.340000, 0.261800,-0.113800,-0.008420},
								 { 6.1510,  0.4524, 83.0000, 18.3100,  2.7100, 4.361000, 0.255900,-0.114500,-0.008447},
								 { 6.5700,  0.4322, 84.7600, 15.5300,  2.7790, 4.349000, 0.240000,-0.116600,-0.008550},
								 { 5.7380,  0.4492, 84.6100, 14.1800,  3.1010, 4.362000, 0.232700,-0.117400,-0.008588},
								 { 5.0310,  0.4707, 85.5800, 16.5500,  3.2110, 4.375000, 0.225300,-0.118500,-0.008648},
								 { 4.3200,  0.4947, 76.1400, 10.8500,  5.4410, 4.362000, 0.206900,-0.121400,-0.008815},
								 { 4.6520,  0.4571, 80.7300, 22.0000,  4.9520, 4.346000, 0.185700,-0.124900,-0.009021},
								 { 3.1140,  0.5236, 76.6700,  7.6200,  6.3850, 4.355000, 0.180000,-0.125500,-0.009045},
								 { 3.1140,  0.5236, 76.6700,  7.6200,  7.5020, 4.389000, 0.180600,-0.125300,-0.009028},
								 { 3.1140,  0.5236, 76.6700,  7.6200,  8.5140, 4.407000, 0.175900,-0.125800,-0.009054},
								 { 5.7460,  0.4662, 79.2400,  1.1850,  7.9930, 4.419000, 0.169400,-0.126700,-0.009094},
								 { 2.7920,  0.6346,106.1000,  0.2986,  2.3310, 4.412000, 0.154500,-0.128900,-0.009202},
								 { 4.6670,  0.5095,124.3000,  2.1020,  1.6670, 4.419000, 0.144800,-0.130300,-0.009269},
								 { 2.4400,  0.6346,105.0000,  0.8300,  2.8510, 4.436000, 0.144300,-0.129900,-0.009229},
								 { 1.4910,  0.7118,120.6000,  1.1010,  1.8770, 4.478000, 0.160800,-0.126200,-0.008983},
								 {11.7200,  0.3826,102.8000,  9.2310,  4.3710, 4.489000, 0.151700,-0.127800,-0.009078},
								 { 7.1260,  0.4804,119.3000,  5.7840,  2.4540, 4.514000, 0.155100,-0.126800,-0.009005},
								 {11.6100,  0.3955,146.7000,  7.0310,  1.4230, 4.533000, 0.156800,-0.126100,-0.008945},
								 {10.9900,  0.4100,163.9000,  7.1000,  1.0520, 4.548000, 0.157200,-0.125600,-0.008901},
								 { 9.2410,  0.4275,163.1000,  7.9540,  1.1020, 4.553000, 0.154400,-0.125500,-0.008883},
								 { 9.2760,  0.4180,157.1000,  8.0380,  1.2900, 4.548000, 0.148500,-0.125900,-0.008889},
								 { 3.9990,  0.6152, 97.6000,  1.2970,  5.7920, 4.489000, 0.112800,-0.130900,-0.009107},
								 { 4.3060,  0.5658, 97.9900,  5.5140,  5.7540, 4.402000, 0.066560,-0.137500,-0.009421},
								 { 3.6150,  0.6197, 86.2600,  0.3330,  8.6890, 4.292000, 0.010120,-0.145900,-0.009835},
								 { 5.8000,  0.4900,147.2000,  6.9030,  1.2890, 4.187000,-0.045390,-0.154200,-0.010250},
								 { 5.6000,  0.4900,130.0000, 10.0000,  2.8440, 4.577000, 0.130000,-0.128500,-0.009067},
								 { 3.5500,  0.6068,124.7000,  1.1120,  3.1190, 4.583000, 0.125300,-0.129100,-0.009084},
								 { 3.6000,  0.6200,105.8000,  0.1692,  6.0260, 4.580000, 0.117400,-0.130100,-0.009129},
								 { 5.4000,  0.5300,103.1000,  3.9310,  7.7670, 4.581000, 0.111000,-0.130900,-0.009161},
								 { 3.9700,  0.6459,131.8000,  0.2233,  2.7230, 4.582000, 0.104600,-0.131700,-0.009193},
								 { 3.6500,  0.6400,126.8000,  0.6834,  3.4110, 4.600000, 0.105200,-0.131500,-0.009178},
								 { 3.1180,  0.6519,164.9000,  1.2080,  1.5100, 4.614000, 0.104300,-0.131500,-0.009175},
								 { 2.0310,  0.7181,153.1000,  1.3620,  1.9580, 4.619000, 0.097690,-0.132500,-0.009231},
								 {14.4000,  0.3923,152.5000,  8.3540,  2.5970, 4.671000, 0.113600,-0.129800,-0.009078},
								 {10.9900,  0.4599,138.4000,  4.8110,  3.7260, 4.706000, 0.120600,-0.128700,-0.009009},
								 {16.6000,  0.3773,224.1000,  6.2800,  0.9121, 4.732000, 0.124400,-0.128000,-0.008968},
								 {10.5400,  0.4533,159.3000,  4.8320,  2.5290, 4.722000, 0.115600,-0.129200,-0.009030},
								 {10.3300,  0.4502,162.0000,  5.1320,  2.4440, 4.710000, 0.106000,-0.130500,-0.009100},
								 {10.1500,  0.4471,165.6000,  5.3780,  2.3280, 4.698000, 0.096470,-0.131900,-0.009169},
								 { 9.9760,  0.4439,168.0000,  5.7210,  2.2580, 4.681000, 0.085360,-0.133500,-0.009252},
								 { 9.8040,  0.4408,176.2000,  5.6750,  1.9970, 4.676000, 0.078190,-0.134500,-0.009302},
								 {14.2200,  0.3630,228.4000,  7.0240,  1.0160, 4.663000, 0.068670,-0.135800,-0.009373},
								 { 9.9520,  0.4318,233.5000,  5.0650,  0.9244, 4.676000, 0.068610,-0.135700,-0.009363},
								 { 9.2720,  0.4345,210.0000,  4.9110,  1.2580, 4.649000, 0.053620,-0.137900,-0.009480},
								 {10.1300,  0.4146,225.7000,  5.5250,  1.0550, 4.634000, 0.043350,-0.139400,-0.009558},
								 { 8.9490,  0.4304,213.3000,  5.0710,  1.2210, 4.603000, 0.026790,-0.141800,-0.009690},
								 {11.9400,  0.3783,247.2000,  6.6550,  0.8490, 4.584000, 0.014940,-0.143600,-0.009783},
								 { 8.4720,  0.4405,195.5000,  4.0510,  1.6040, 4.576000, 0.007043,-0.144700,-0.009841},
								 { 8.3010,  0.4399,203.7000,  3.6670,  1.4590, 4.571000, 0.000705,-0.145600,-0.009886},
								 { 6.5670,  0.4858,193.0000,  2.6500,  1.6600, 4.566000,-0.005626,-0.146400,-0.009930},
								 { 5.9510,  0.5016,196.1000,  2.6620,  1.5890, 4.561000,-0.011970,-0.147300,-0.009975},
								 { 7.4950,  0.4523,251.4000,  3.4330,  0.8619, 4.572000,-0.012000,-0.147200,-0.009965},
								 { 6.3350,  0.4825,255.1000,  2.8340,  0.8228, 4.569000,-0.017550,-0.148000,-0.010000},
								 { 4.3140,  0.5558,214.8000,  2.3540,  1.2630, 4.573000,-0.019920,-0.148200,-0.010010},
								 { 4.0200,  0.5681,219.9000,  2.4020,  1.1910, 4.570000,-0.025470,-0.149000,-0.010050},
								 { 3.8360,  0.5765,210.2000,  2.7420,  1.3050, 4.528000,-0.046130,-0.152100,-0.010220},
								 { 4.6800,  0.5247,244.7000,  2.7490,  0.8962, 4.494000,-0.063700,-0.154800,-0.010370},
								 { 3.2230,  0.5883,232.7000,  2.9540,  1.0500, 4.564000,-0.027000,-0.147100,-0.009852},
								 { 2.8920,  0.6204,208.6000,  2.4150,  1.4160, 4.546000,-0.049630,-0.152300,-0.010220},
								 { 4.7280,  0.5522,217.0000,  3.0910,  1.3860, 4.594000,-0.033390,-0.149600,-0.010060},
								 { 6.1800,  0.5200,170.0000,  4.0000,  3.2240, 4.608000,-0.028860,-0.148500,-0.009990},
								 { 9.0000,  0.4700,198.0000,  3.8000,  2.0320, 4.624000,-0.026390,-0.148100,-0.009971},
								 { 2.3240,  0.6997,216.0000,  1.5990,  1.3990, 4.636000,-0.024220,-0.147700,-0.009939},
								 { 1.9610,  0.7286,223.0000,  1.6210,  1.2960, 4.648000,-0.021720,-0.147100,-0.009903},
								 { 1.7500,  0.7427,350.1001,  0.9789,  0.5507, 4.662000,-0.119200,-0.175200,-0.011960},
								 {10.3100,  0.4613,261.2000,  4.7380,  0.9899, 4.690000,-0.009867,-0.144900,-0.009771},
								 { 7.9620,  0.5190,235.7000,  4.3470,  1.3130, 4.715000,-0.002113,-0.143500,-0.009689},
								 { 6.2270,  0.5645,231.9000,  3.9610,  1.3790, 4.729000, 0.001392,-0.142800,-0.009644},
								 { 5.2460,  0.5947,228.6000,  4.0270,  1.4320, 4.729000,-0.000598,-0.143000,-0.009647},
								 { 5.4080,  0.5811,235.7000,  3.9610,  1.3580, 4.738000, 0.001075,-0.142500,-0.009618},
								 { 5.2180,  0.5828,245.0000,  3.8380,  1.2500, 4.751000, 0.004244,-0.141900,-0.009576}};

	return res[z-1][index-1];
}

double Kinematic::HeliumB(int index, int z) {
	double res[92][5] = {{ 0.39, 0.63, 4.17, 85.55, 19.55},
								{ 0.58, 0.59, 6.3, 130., 44.07},
								{ 15.23, 0.1076, 10.14, 1232., 31.24},
								{ 2.2060,  0.5100, 15.3200,  0.2500,  8.9950},
								{ 3.6910,  0.4128, 18.4800, 50.7200,  9.0000},
								{ 3.47,0.4485, 22.37, 36.41, 7.993},
								{ 2.0, 0.548, 29.82, 18.11, 4.37},
								{ 2.717, 0.4858, 32.88, 25.88, 4.336},
								{ 2.616, 0.4708, 41.2, 28.07, 2.458},
								{ 2.303, 0.4861, 37.01, 37.96, 5.092},
								{ 13.03, 0.2685, 35.65, 44.18, 9.175},
								{ 4.3000,  0.4700, 34.3000,  3.3000, 12.7400},
								{ 2.5000,  0.6250, 45.7000,  0.1000,  4.3590},
								{ 2.1000,  0.6500, 49.3400,  1.7880,  4.1330},
								{ 1.7290,  0.6562, 53.4100,  2.4050,  4.8450},
								{ 3.116,   0.5988, 53.71,    5.632,    4.536},
								{ 3.365,   0.571,  63.67,    6.182,    2.969},
								{ 2.291,   0.6284, 73.88,    4.478,    2.066},
								{ 16.6,    0.3095, 99.1,     10.98,    1.092},
								{ 6.2970,  0.4622, 65.3900, 10.1400,  5.0360},
								{ 5.3070,  0.4918, 61.7400, 12.4000,  6.6650},
								{ 4.7100,  0.5087, 65.2800,  8.8060,  5.9480},
								{ 6.1510,  0.4524, 83.0000, 18.3100,  2.7100},
								{ 6.5700,  0.4322, 84.7600, 15.5300,  2.7790},
								{ 5.7380,  0.4492, 84.6100, 14.1800,  3.1010},
								{ 5.0310,  0.4707, 85.5800, 16.5500,  3.2110},
								{ 4.3200,  0.4947, 76.1400, 10.8500,  5.4410},
								{ 4.6520,  0.4571, 80.7300, 22.0000,  4.9520},
								{ 3.1140,  0.5236, 76.6700,  7.6200,  6.3850},
								{ 3.1140,  0.5236, 76.6700,  7.6200,  7.5020},
								{ 3.1140,  0.5236, 76.6700,  7.6200,  8.5140},
								{ 5.7460,  0.4662, 79.2400,  1.1850,  7.9930},
								{ 2.7920,  0.6346,106.1000,  0.2986,  2.3310},
								{ 4.6670,  0.5095,124.3000,  2.1020,  1.6670},
								{ 1.65,    0.7,    148.1,      1.47,  0.9686},
								{ 1.413,   0.7377, 147.9,     1.466,  1.0160},
								{11.7200,  0.3826,102.8000,  9.0000,  4.3710},
								{ 7.1260,  0.4804,119.3000,  5.7840,  2.4540},
								{11.6100,  0.3955,146.7000,  7.0310,  1.4230},
								{10.9900,  0.4100,163.9000,  7.1000,  1.0520},
								{ 9.2410,  0.4275,163.1000,  7.9540,  1.1020},
								{ 9.2760,  0.4180,157.1000,  8.0380,  1.2900},
								{ 3.9990,  0.6152, 97.6000,  1.2970,  5.7920},
								{ 4.3020,  0.5658, 97.9900,  5.5140,  5.7540},
								{ 3.6150,  0.6197, 86.2600,  0.3330,  8.6890},
								{ 5.8000,  0.4900,147.2000,  6.9030,  1.2890},
								{ 5.6000,  0.4900,130.0000, 10.0000,  2.8440},
								{ 3.5500,  0.6068,124.7000,  1.1120,  3.1190},
								{ 3.6000,  0.6200,105.8000,  0.1692,  6.0260},
								{ 5.4000,  0.5300,103.1000,  3.9310,  7.7670},
								{ 3.9700,  0.6459,131.8000,  0.2233,  2.7230},
								{ 3.6500,  0.6400,126.8000,  0.6834,  3.4110},
								{ 4.13,    0.6177, 152.0,    2.516,   1.9380},
								{ 3.949,   0.6209, 200.5,    1.878,   0.9126},
								{ 25.94,   0.3139, 335.1,    2.946,   0.3347},
								{710.9900,  0.4599,138.4000, 4.8110,  3.7260},
								{16.6000,  0.3773,224.1000,  6.2800,  0.9121},
								{10.5400,  0.4533,159.3000,  4.8320,  2.5290},
								{10.3300,  0.4502,162.0000,  5.1320,  2.4440},
								{10.1500,  0.4471,165.6000,  5.3780,  2.3280},
								{ 9.9760,  0.4439,168.0000,  5.7210,  2.2580},
								{ 9.8040,  0.4408,176.2000,  5.6750,  1.9970},
								{14.2200,  0.3630,228.4000,  7.0240,  1.0160},
								{ 9.9520,  0.4318,233.5000,  5.0650,  0.9244},
								{ 9.2720,  0.4345,210.0000,  4.9110,  1.2580},
								{10.1300,  0.4146,225.7000,  5.5250,  1.0550},
								{ 8.9490,  0.4304,213.3000,  5.0710,  1.2210},
								{11.9400,  0.3783,247.2000,  6.6550,  0.8490},
								{ 8.4720,  0.4405,195.5000,  4.0510,  1.6040},
								{ 8.3010,  0.4399,203.7000,  3.6670,  1.4590},
								{ 6.5670,  0.4858,193.0000,  2.6500,  1.6600},
								{ 5.9510,  0.5016,196.1000,  2.6620,  1.5890},
								{ 7.4950,  0.4523,251.4000,  3.4330,  0.8619},
								{ 6.3350,  0.4825,255.1000,  2.8340,  0.8228},
								{ 4.3140,  0.5558,214.8000,  2.3540,  1.2630},
								{ 4.0200,  0.5681,219.9000,  2.4020,  1.1910},
								{ 3.8360,  0.5765,210.2000,  2.7420,  1.3050},
								{ 4.6800,  0.5247,244.7000,  2.7490,  0.9862},
								{ 3.2230,  0.5883,232.7000,  2.9540,  1.0500},
								{ 8.15,    0.4745, 269.2,    2.392,   0.7467},
								{ 4.7280,  0.5522,217.0000,  3.0910,  1.3860},
								{ 6.1800,  0.5200,170.0000,  4.0000,  3.2240},
								{ 9.0000,  0.4700,198.0000,  3.8000,  2.0320},
								{ 2.3240,  0.6997,216.0000,  1.5990,  1.3990},
								{ 1.9610,  0.7286,223.0000,  1.6210,  1.2960},
								{ 4.822,   0.605, 418.3,     0.8335,  0.3865},
								{10.3100,  0.4613,261.2000,  4.7380,  0.9899},
								{ 7.9620,  0.5190,235.7000,  4.3470,  1.3130},
								{ 6.2270,  0.5645,231.9000,  3.9610,  1.3790},
								{ 5.2460,  0.5947,228.6000,  4.0270,  1.4320},
								{ 5.4080,  0.5811,235.7000,  3.9610,  1.3580},
								{ 5.2180,  0.5828,245.0000,  3.8380,  1.2500}};

	return res[z-1][index-1];
}

double Kinematic::Momentum(double Ekin, double mass) {
	return momentum(Ekin, mass);
}

/************************************
  This file contains the code for doing kinematics
  boosts rotations....
 ******************************************/

double gen_theta();

/**********************************************
  v = dx/dt
  dp/dt = qv x B(x)

  v is velocity in units of c
  q is the charge in units of e
  B is magnetic field in Tesla
  p is the momentum in MeV/c
  t is in units of (s) 
  masses are given in MeV/c2

 ***********************************************/
double dot(double* x,double* y) {
	/*     double *x,*y; */
	int i;
	double retval=0.0;

	for(i=0;i<3;i++) {
		retval += (*x)*(*y);
		x++;
		y++;
	}
	return(retval);
}

double norm(double* x) {
	/*     double *x; */
	/************************************************************
	  Calculate the norm of a three vector
	 **************************************************************/
	double retval;

	retval = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	retval = sqrt(retval);
	return(retval);
}

void cross_product(double* a, double* b, double* c) {
	/********     double *a,*b,*c; */
	/**********************************************************
	  calculate a = b x c
	 **********************************************************/
	a[0] = b[1]*c[2]-b[2]*c[1];
	a[1] = b[2]*c[0]-b[0]*c[2];
	a[2] = b[0]*c[1]-b[1]*c[0];
}

void scaler_mul(double a,double* v) {
	/*     double a,*v; */
	/**************************************************************
	  ( a*v[] )=> v[]
	  Replace the vector v with a*v
	 ***************************************************************/
	*v++ *= a;
	*v++ *= a;
	*v++ *= a;
}


double momentum(double ke,double m) {
	/*     double ke,m; */
	/*************************************************************************
	  E = ke + m
	  and 
	  E*E = p*p + m*m
	  p = sqrt((ke+m)^2 - m*m)
	  = sqrt(ke*ke + 2*ke*m)
	  ke and m must be in units of MeV
	  p is returned in units of MeV/c
	 ***************************************************************************/
	double p;
	p  = sqrt(ke*ke + ((double)2.0)*ke*m);
	return(p);
}

double speed(double p,double m) {
	/*     double p,m; */
	/********************************************************************************
	  return the speed of a particle
	  given its momentum and mass
	 ******************************************************************************/
	double retval;
	retval = (p/m)/sqrt((double)1.0+ (p*p)/(m*m));
#ifdef orbits_debug_results
	std::cout<<"speed("<<p<<", "<<m<<"): "<<retval<<" = ("<<p<<"/"<<m<<")/sqrt(1.0+"<<p*p<<"/"<<m*m<<")"<<std::endl;
#endif
	return(retval);
}

void velocity(double* v,double* p,double m) {
	/*     double *v,*p,m; */
	/************************************************************************
	  take p[] and m
	  and calculate v[]
	  p is units of MeV/c
	  m is units of MeV/c2
	  returns v in units of c
	 *************************************************************************/
	double mag_p,mag_v;
	int i;

	mag_p = norm(p);
	mag_v = speed(mag_p,m);
	for (i=0;i<3;i++) {
		v[i] = (mag_v*p[i])/mag_p;
	}
}

double rel_gamma(double v) {
	/*     double v; */
	/**********************************************************************************
	  v is the speed in units where c=1
	 *********************************************************************************/
	double g=1.0/sqrt(1.0-v*v);
	return(g);
}

void boost(double* v,double* b) {
	/*     double *v,*b; */
	/********************************************************
	  v is replaced by boosted velocity
	 ********************************************************/
	double vp,vt,vp_prime,vt_prime,dummy,dvb,bmag;
	int i;


	bmag=norm(b);
	if(bmag > 0.0 ) {
		dvb =dot(v,b);
		vp = dvb/bmag;
		dummy= norm(v);
		vt = sqrt(dummy*dummy-vp*vp);
		vp_prime=(vp-bmag)/(1.0-dvb);
		vt_prime=vt/(rel_gamma(bmag)*(1.0-dvb));
		for(i=0;i<3;i++) {
			/* this constructs  perpendicular vector */
			v[i] -= b[i]*vp/bmag;
		}
		dummy = norm(v);
		for(i=0;i<3;i++) {
			if(dummy >0.0) {
				v[i] *= vt_prime/dummy;
			} else {
				v[i] = 0.0;
			}
			v[i] += b[i]*vp_prime/bmag;
		}
	}
}


double BetaCm(double Ebeam,Element** particle) {
	/*     double Ebeam;
			 ISOTOPEPTR  *particle;*/

	/************************************
	  particle[0] is the target
	  particle[1] is the beam
	  Ebeam is beam energy in MeV
	  returns the Beta of the COMS 
	 ***********************************/
	double p1,bcm,m1,m0;

	m1=particle[1]->Mass();
	m0=particle[0]->Mass();
	p1 = momentum(Ebeam,m1);
	bcm = p1/(Ebeam+m1+m0);
#ifdef debug_orbits_results
	std::cout<<"BetaCm("<<Ebeam<<", ...): "<<bcm<<" = "<<p1<<"/("<<Ebeam<<"+"<<m1<<"+"<<m0<<")"<<std::endl;
#endif
	return(bcm);

}


double EnergyCm(double beam_energy, Element** particle) {
	/*     double Ebeam;
			 ISOTOPEPTR  *particle;  */

	/************************************
	  particle[0] is the target
	  particle[1] is the beam
	  beam_energy is beam energy in MeV
	  Returns Total Energy in COMS
	 ***********************************/

	double w,m1,m0;
	m1=particle[1]->Mass();
	m0=particle[0]->Mass();
	w = m1*m1+m0*m0+2.0*m0*(m1+beam_energy);
	w = sqrt(w);
#ifdef debug_orbits_results
	std::cout<<"EnergyCm("<<beam_energy<<", ...): "<<w<<" = sqrt("<<m1*m1<<"+"<<m0*m0<<"+ 2.0*"<<m0<<"*("<<m1<<"+"<<beam_energy<<"))"<<std::endl;
#endif
	return(w);
}

double MomentumCm(double ecm,double ex2,double ex3,Element** particle) {
	/*     double ecm,ex2,ex3;
			 ISOTOPEPTR *particle; */
	/****************************
	  particle[2] is the detected ejectile
	  particle[3] is the recoil particle
	  ex2 the excitation of particle[2]
	  ex3 the excitation of particle[3]
	  ecm is the total energy in the CM = m0+m1+Ebeam
	  returns the COMS 3-momentum magnitiude of the scattered particles

	  error conditions:

	  will return a -1.0 if ecm<threshold energy
	 ***************************************/

	double a,b,c,q,m2,m3,d;

	m2=particle[2]->Mass()+ex2;
	m3=particle[3]->Mass()+ex3;
	a = ecm*ecm;
	b = a;
	a *= a;
	c = m2*m2+m3*m3;
	d = m2*m2-m3*m3;
	q = 0.25*(a+d*d-2.0*b*c)/b;
	if (q<0.0) return((double)(-1.0));
	q  = sqrt(q);
#ifdef debug_orbits_results
	std::cout<<"MomentumCm("<<ecm<<", "<<ex2<<", "<<ex3<<", ...): "<<q<<" = sqrt(0.25*("<<a<<"+"<<d*d<<"- 2.0*"<<b<<"*"<<c<<")/"<<b<<")"<<std::endl;
#endif
	return(q);
}

double mom_from_speed(double v,double m) {
	/*     double v,m; */
	double p;
	p= v*m/sqrt(1.0-v*v);
	return(p);
}

double theta_cm(double tlab,double bcm,double vlab) {
	/*     double tlab,bcm,vlab; */

	double t1;

	t1 = vlab*sin(tlab)/(rel_gamma(bcm)*(vlab*cos(tlab)-bcm));
	t1 = atan(t1);
	if ( t1<0.0) t1 +=TMath::Pi();
	return(t1);
}

double theta_lab_max(double ebeam,double ex2,double ex3,Element** particle) {
	/*     double ebeam,ex2,ex3;
			 ISOTOPEPTR *particle;  */
	double ecm,b3[3],v3[3],u,p,m,ucm,t_max_cm,temp;
	double ret_val;

	ecm=EnergyCm(ebeam,particle);
	p = MomentumCm(ecm,ex2,ex3,particle);
	m=particle[2]->Mass()+ex2;
	u = speed(p,m);
	ucm = BetaCm(ebeam,particle);
	if (u>ucm) {
		return((double)TMath::Pi());
	}
	else {
		t_max_cm=acos(-1.0*u/ucm);
		v3[0] = cos(t_max_cm)*u;
		v3[1] = sin(t_max_cm)*u;
		v3[2] = 0.0;
		b3[0] = -1.0*ucm;
		b3[1] = 0.0;
		b3[2] = 0.0;
		boost(v3,b3);
		temp = v3[0]/norm(v3);
		ret_val = acos(temp);
		return(ret_val);
	}
}

int scatter(double* p2,double* p3,
		double Ebeam,
		double ex2,double ex3,
		double* cm_unit,Element** particle) {
	/*     double *p2,*p3,Ebeam,*cm_unit,ex2,ex3;
			 ISOTOPEPTR *particle;  */

	/******************************************
	  Routine takes particle[1] at Lab Energy Ebeam
	  and scatters it off Particle[0] into  paarticles[2] and [3]
	  at center  of mass angle cm_angle.
	  Where ex2 and ex3 are the respective out-going.
	 ************************************************/

	double v,b[3],v2[3],v3[3],p,m3,m2;
	double dummy;
	int i;

	m2 = particle[2]->Mass()+ex2;
	m3 = particle[3]->Mass()+ex3;
	dummy = EnergyCm(Ebeam,particle);
	if(dummy<0.0) return(0);
	p = MomentumCm(dummy,ex2,ex3,particle);
	b[2] = -1.0*BetaCm(Ebeam,particle);
	b[1] = 0.0;
	b[0] = 0.0;
	v = speed(p,m2);
	for(i=0;i<3;i++) {
		v2[i] = cm_unit[i]*v;
	}
	v = -1.0*speed(p,m3);
	for(i=0;i<3;i++) {
		v3[i] = cm_unit[i]*v;
	}
	boost(v3,b);
	boost(v2,b);
#ifdef debug_orbits_results
	std::cout<<"scatter(): "<<m2<<" "<<m3<<" "<<dummy<<" "<<p<<", ("<<v2[0]<<", "<<v2[1]<<", "<<v2[2]<<"), ("<<v3[0]<<", "<<v3[1]<<", "<<v3[2]<<")"<<std::endl;
#endif
	v_to_p(v3,p3,m3);
	v_to_p(v2,p2,m2);
	return(1);
}

void v_to_p(double* v,double* p,double m) {
	/*     double *v,*p,m;  */
	double mag_v,mag_p;
	int i;
	mag_v = norm(v);
	mag_p = mom_from_speed(mag_v,m);
	if(mag_v == 0.0) {
		p[0]=p[1]=p[2]=0.0;
	}
	else {
		for(i=0;i<3;i++) {
			*p = (*v)*mag_p/mag_v;
			p++;
			v++;
		}
	}
}

double energy(double p,double m) {
	/*     double p,m;  */
	double ke;
	ke = sqrt(p*p+m*m) - m;
	return(ke);
}

double gen_event(double* p2,double* p3,double Ebeam,double ex2,double ex3,Element** particle,double dangle) {
	/*     double *p2,*p3,Ebeam,ex2,ex3,dangle;
			 ISOTOPEPTR *particle;  */
	/*********************************************
	  this routine generates a distributed random
	  event at the origin using a cross section from gen_theta()
	  and boosts it subject to the beam direction
	  give by beam_dir(b) where dangle is the maximum 
	  divergence angle of the beam
	 ******************************************************/

	double v,b[3],v2[3],v3[3],p,m3,m2,cm_unit[3],bmag;
	double dummy;
    int i;
	double theta,phi;

	m2 = particle[2]->Mass()+ex2;
	m3 = particle[3]->Mass()+ex3;
	dummy = EnergyCm(Ebeam,particle);
	if(dummy<0.0) return(-999.0);
	p = MomentumCm(dummy,ex2,ex3,particle);
	bmag = -1.0*BetaCm(Ebeam,particle);
	beam_dir(b,dangle);
	scaler_mul(bmag,b);
	theta = gen_theta();
	phi = 2.0*TMath::Pi()*ran1();
	cm_unit[2] = cos(theta);
	cm_unit[1] = -1.0*sin(theta)*sin(phi);
	cm_unit[0] = sin(theta)*cos(phi);
	v = speed(p,m2);
	for(i=0;i<3;i++) {
		v2[i] = cm_unit[i]*v;
	}
	v = -1.0*speed(p,m3);
	for(i=0;i<3;i++) {
		v3[i] = cm_unit[i]*v;
	}
	boost(v3,b);
	boost(v2,b);
	v_to_p(v3,p3,m3);
	v_to_p(v2,p2,m2);
	return(theta);
}

void beam_dir(double* beam,double thetamax) {
	/***********************************************************
	  generate a unform distribution in solid angle in a cone thetamax
	  about the z-axis.
	 ***********************************************************/
	/*     double *beam,thetamax;  */
	double theta,phi;
	phi = 2.0*TMath::Pi()*ran1();
	theta = acos(1.0 -(1.0-cos(thetamax))*ran1());
	beam[2] = cos(theta);
	beam[1] = (-1.0*sin(phi))*sin(theta);
	beam[0] = (cos(phi)*sin(theta));
}

void gen_offset(double* offset,double xsize,double ysize) {
	/*     double *offset,xsize,ysize;  */
	/**************************************************************
	  generate a uniform distribution over an ellipse of xsize,ysize
	 **************************************************************/
	double x,y,r;

	offset[2] = 0.0;
	r=2;
	while(r>1.0) {
		x = ysize*(2.0*ran1()-1.0);
		y = xsize*(2.0*ran1()-1.0);
		r = x*x/(xsize*xsize)+y*y/(ysize*ysize);
	}
	offset[1] = y;
	offset[0] = x;

}

double gasdev() {
	/***********************************************************
	  this routine returns a double randomly distributed in a guassian
	  of unit norm and zero mean.
	 ************************************************************/
	static int iset = 0;
	static double gset;
	double fac,r,v1,v2;
	//double ran1();

	if (iset == 0) {
		do {
			v1 = 2.0*ran1()-1.0;
			v2 = 2.0*ran1()-1.0;
			r = v1*v1+v2*v2;
		} while (r >=1.0);
		fac = sqrt(-2.0*log(r)/r);
		gset = v1*fac;
		iset = 1;
		return(v2*fac);
	} else {
		iset = 0;
		return(gset);
	}
}

double ran1() {
	/*  int rand(); */
	double ret_val;
#ifndef RAND_MAX
	fprint("RAND_MAX not defined")
#endif
		unsigned long rand_max = RAND_MAX;
	static int itest=0;
	int i,j;
	static double v[97]; 
	if(itest==0) {
		itest=1;
		for(i=0;i<97;i++) {
			v[i] =((double)rand())/((double)rand_max);
		}
	}
	j = rand()*97/rand_max;
	ret_val = v[j];
	v[j]=((double)rand())/((double)rand_max);
	return(ret_val);
} 

double lordev() {
	/***********************************************************

	  Returns a lorentzian deviate of unit width and zero mean.

	 ************************************************************/

	double retval;

	retval = tan(TMath::Pi()*ran1());
	return(retval);
}

double expdev() {
	/*********************************************************
	  returns an exponential deviate of index 1
	  ie exp(-y)
	 **********************************************************/
	double retval;
	retval = -log(ran1());
	return(retval);
}	

double gen_theta() {
	double theta;

	/************************

	  iso_tropic
	  theta = acos(2.0*ran1()-1.0);

	  do/dw = 1/sin(theta)
	  theta=ran1()*TMath::Pi();
	 ************************/
	theta=ran1()*TMath::Pi()/2.0;		/**** 1/sin(theta) limited to forward 90 deg****/
	return(theta);
}












