#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <fstream>
#include <valarray>
#include <cmath>
#include <iomanip>
#include <vector>
//#include <omp.h>
#include <sys/time.h>
#include <iostream>
//#include "config.h"




#include <TROOT.h>

#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TString.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TSpline.h>

#include <TComplex.h>
#include "TTimer.h"

using namespace std;


// ##### MODEL CONFIGURATION
constexpr double    dt = 0.5e-1;
//constexpr int     PARN = 3*2;// 500;
constexpr int     MAX_PAR_nPart = 100;//20;//3;//5;//2;//7;//2;//2;//3;//2; // N particles
int     nPart = 3; // N particles set from outside
//constexpr int     PAR_nDim = 3; // nDimensions
//constexpr double    m = 1.;
int nVar = 6*nPart;

double masses[MAX_PAR_nPart];
int charges[MAX_PAR_nPart];


// for termination of steps:
bool global_flag_check_distances = false;
double prev_dr[MAX_PAR_nPart][MAX_PAR_nPart];
double MAX_DIST_BETWEEN_PARTICLES = 4;//2;//4; // fm
bool global_terminate_flag = false;

//constexpr double PARka = 2.0;
//constexpr double PARka_over_N = PARka/(double)PARN;
//constexpr double PARw = M_PI;
constexpr auto num_steps = 100;


bool PRINT_DEBUG = false;
bool TXT_OUTPUT_TO_FILE = false;


double minPotRange = 0.1;//35;//5;//0.3;
double maxPotRange = 10.0;



// ##### POTENTIALS:
const int N_POT_POINTS = 1000;
TF1 *fPotReidNoCoulumb;
TF1 *fPotReidCoulumbAttract;
TF1 *fPotReidCoulumbRepulsive;
TF1 *fPotCoulumbOnly;

TSpline5 *spline5_Av18 = 0;
TSpline5 *spline5_Av18_ATTRACTIVE = 0;
TSpline5 *spline5_Av18_REPULSIVE  = 0;


TSpline5 *spline5_FORCE_Av18_ATTRACTIVE = 0;
TSpline5 *spline5_FORCE_Av18_REPULSIVE  = 0;

TGraphErrors *grF_attr;
TGraphErrors *grF_rep;
TGraphErrors *grF_Av18;
TGraphErrors *grF_Av18_ATTRACTIVE;
TGraphErrors *grF_Av18_REPULSIVE;


ofstream *outTxtFile = 0x0;



// ##### FOR ANIMATION DRAWING:
bool DO_ANIMATION = false;//true;

//TEllipse *elPart[MAX_PAR_nPart][num_steps];
//int timeStepCounter = 0;
TEllipse *el_IS_Part[MAX_PAR_nPart];
TEllipse *elAnimation[MAX_PAR_nPart];
TLine *lineAnimation[MAX_PAR_nPart*MAX_PAR_nPart];

TPolyMarker3D *pmAnimation[MAX_PAR_nPart];
TPolyLine3D *plAnimation[MAX_PAR_nPart*MAX_PAR_nPart];

double systemRadius = 15;//8;
double fVisNucleusRadiusNucleus = 0.05; //visual size nucleus
double kNucleonSize = 0.8;
double rVisNucleon = fVisNucleusRadiusNucleus * kNucleonSize/systemRadius;


void animRunner();


// other model declarations:
using StateVec   = std::valarray<double>;
using DerivFunc  = StateVec (*) (const StateVec &);
using OutputFunc = void (*) (double , const StateVec &, std::ofstream &);
using Integrator = void (*) (StateVec &, DerivFunc, double);

//StateVec Kuramoto (const StateVec &x);
StateVec Coulomb (const StateVec &x);
void euler_integrator (StateVec &, DerivFunc, double);
void runge_kutta4_integrator (StateVec &, DerivFunc, double);
void integrate (Integrator, DerivFunc, OutputFunc, 
                /*std::ofstream &*/ std::ofstream *output, StateVec, int, double);

void print_coordinates (double time, const StateVec &y, 
                        std::ofstream &output);
double get_wall_time();
double get_cpu_time();
void display_timing(double , double );

void init();
int runRK(StateVec &outsideState); // , bool standaloneUse , int nPart);
void runRKstandalone() ;

StateVec *yAllSteps[num_steps+1];
StateVec yFinalState;



#endif
