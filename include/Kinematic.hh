#ifndef KINEMATIC
#define KINEMATIC

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>

#include "TMath.h"
#include "TSpline.h"

#include "Element.hh"
#include "Material.hh"
#include "Isotopes.hh"

//#include "Utilities.hh"


#ifndef PI
#define PI                       (TMath::Pi())
#endif

class Kinematic {
public:
  Kinematic();
  Kinematic(Element* pro, Material* tar);
  Kinematic(Element* pro, Material* tar, double target_thick);
  ~Kinematic();

  int orbits(Element* projectile, Element* target, Element* recoil, Element* ejectile, double energy, double cm_start, double cm_stop, double cm_step, double ex_start, double ex_stop, double ex_step, int code, double**& results_en, double**& results_ang);
  int orbits(Element* projectile, Element* target, Element* recoil, Element* ejectile, double energy, double cm_start, double cm_stop, double cm_step, double ex, int code, double*& results_en, double*& results_ang);
  // included by S.Klupp
  int orbits(Element* projectile, Element* target, Element* recoil, Element* ejectile, double energy, double cm, double ex, int code, double& results_en, double& results_ang);
  int orbits(Element* projectile, Element* target, Element* recoil, Element* ejectile, double energy, double cm_start, double cm_stop, double cm_step, double ex, int code, std::vector<double>& results_en, std::vector<double>& results_ang);

  double CMAngle(double LabAngle, Element* projectile, Element* target, Element* recoil, Element* ejectile, int energies, double* energy, double ex, int code, double limit = 1e-3);
  void SetThickness(double target_thick);
  void Projectile(Element*);
  void Target(Material*);
  int Run(double energy, double*& results);
  double Range(double);
  double EnergyLoss(double);
  double EnergyAfterTarget(double);
  double Energy(double,double,double);
  //double ReverseEnergy(double,double,double,double*&);

  TSpline3* RangeVsEnergy(double,double);
  TSpline3* EnergyVsRange(double,double);
  TSpline3* EnergyVsThickness(double,double);
  TSpline3* RutherfordVsThickness(double,double,bool);
  TSpline3* EnergyLossVsEnergy(double, double step_size = 0.1, double factor = 1.);

  void EnergyLossVsEnergy(double, std::vector<double>&, std::vector<double>&, double step_size = 0.1, double factor = 1.);
  double Momentum(double Ekin, double mass);
private:
  double CompoundRange(Element* projectile, Material* target, double energy, int integration_limit);
  double CompoundStoppingPower(Element* projectile, Material* target, double energy);
  double StoppingPower(Element* projectile, Element* target, double energy, bool gaseous);
  double EnergyLoss(Element* projectile, Material* target, double target_thickness, double energy, int integration_limit);
  double EnergyStraggling(Element* projectile, Material* target, double dE_dx, double dE_dx_after_target, double energy_loss);
  double AngularStraggling(Element* projectile, Material* target, double target_thickness, double energy);
  double ChargeState(Element* projectile, double energy);
  double HydrogenA(int index, int z);
  double HeliumA(int index, int z);
  double HeliumB(int index, int z);
  double ShellCorrection(int z);

  Element* fProjectile;
  Material* fTarget;
  double fTargetThickness;
};

double dot(double* x,double* y);
double norm(double* x);
void cross_product(double* a, double* b, double* c);
void scaler_mul(double a,double* v);
double momentum(double ke,double m);
double speed(double p,double m);
void velocity(double* v,double* p,double m);
double rel_gamma(double v);
void boost(double* v,double* b);
double BetaCm(double Ebeam,Element** particle);
double EnergyCm(double beam_energy, Element** particle);
double MomentumCm(double ecm,double ex2,double ex3,Element** particle);
double mom_from_speed(double v,double m);
double theta_cm(double tlab,double bcm,double vlab);
double theta_lab_max(double ebeam,double ex2,double ex3,Element** particle);
int scatter(double* p2,double* p3,double Ebeam,double ex2,double ex3,double* cm_unit,Element** particle);
void v_to_p(double* v,double* p,double m);
double energy(double p,double m);
double gen_event(double* p2,double* p3,double Ebeam,double ex2,double ex3,Element** particle,double dangle);
void beam_dir(double* beam,double thetamax);
void gen_offset(double* offset,double xsize,double ysize);
double gasdev();
double ran1();
double lordev();
double expdev();
double gen_theta();

#endif

