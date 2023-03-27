#include "TROOT.h"
#include "TSystem.h"
//#include "functions.h"


using namespace std;

//int run ();

void runAll()
{
    // #### runge-kutta example is taken from here: https://gitlab.com/a.ziaeemehr/cpp_workshop/tree/master/ODE/ode_oop/Kuramoto
    // #### advised here: https://codereview.stackexchange.com/questions/167460/runge-kutta-fourth-order-integration



//  gROOT->ProcessLine( ".L config.h" );
  gROOT->ProcessLine( ".L functions.cpp" );
  gROOT->ProcessLine( ".L runRK.cpp" );

//  run();

//  return 0;
}
