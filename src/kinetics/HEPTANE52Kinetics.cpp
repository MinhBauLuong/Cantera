/**
 *  @file HEPTANE52Kinetics.cpp 
 *
 * @ingroup chemkinetics
 */

// Author: T. Jaravel

#include "cantera/kinetics/HEPTANE52Kinetics.h"

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
void heptane52_(doublereal* P, doublereal* T, const doublereal* m_y,doublereal* wdot);
}
      
HEPTANE52Kinetics::
HEPTANE52Kinetics(thermo_t* th) : GasKinetics(th) {
//printf("\nWarning: You are using a 19 species reduced scheme for CH4 Air\n");
}

void HEPTANE52Kinetics::get_wdot_reduced(doublereal* wdot)
{
compute_reduced(thermo().pressure(), thermo().temperature(),thermo().massFractions(),wdot);
} 

void HEPTANE52Kinetics::compute_reduced(doublereal P, doublereal T, const doublereal* m_y,doublereal* wdot)
{
int i;
//pressure unit conversion for chemkin
P=P*10.0;
heptane52_(&P,&T,&m_y[0],&wdot[0]);

//mol/mmol conversion cantera is in mmol
for (i=0;i<52;i++) {
	wdot[i]=wdot[i]*1000.0;

}

}

}
