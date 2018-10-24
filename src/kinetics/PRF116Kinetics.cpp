/**
 *  @file PRF116Kinetics.cpp 
 *
 * @ingroup chemkinetics
 */

// Author: T. Jaravel

#include "cantera/kinetics/PRF116Kinetics.h"

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Enhanced3BConc.h"
#include "cantera/kinetics/ThirdBodyMgr.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/thermo/IdealGasPhase.h"
#include <iostream>
 
using namespace std;

namespace Cantera 
{
  //Fortran External Routine
extern "C" 
{ 
void prf116_(doublereal* P, doublereal* T, const doublereal* m_y,doublereal* wdot);
}
      
PRF116Kinetics::
PRF116Kinetics(thermo_t* th) : GasKinetics(th) {
//printf("\nWarning: You are using a 19 species reduced scheme for CH4 Air\n");
}

void PRF116Kinetics::get_wdot_reduced(doublereal* wdot)
{
compute_reduced(thermo().pressure(), thermo().temperature(),thermo().massFractions(),wdot);
} 

void PRF116Kinetics::compute_reduced(doublereal P, doublereal T, const doublereal* m_y,doublereal* wdot)
{
int i;
//pressure unit conversion for chemkin
P=P*10.0;
prf116_(&P,&T,&m_y[0],&wdot[0]);

//mol/mmol conversion cantera is in mmol
for (i=0;i<116;i++) {
	wdot[i]=wdot[i]*1000.0;

}

}

}
