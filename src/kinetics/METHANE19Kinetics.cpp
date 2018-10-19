/**
 *  @file METHANE19Kinetics.cpp 
 *
 * @ingroup chemkinetics
 */

// Author: T. Jaravel

#include "cantera/kinetics/METHANE19Kinetics.h"

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
void methane19_(doublereal* P, doublereal* T, const doublereal* m_y,doublereal* wdot);
}
      
METHANE19Kinetics::
METHANE19Kinetics(thermo_t* th) : GasKinetics(th) {
//printf("\nWarning: You are using a 19 species reduced scheme for CH4 Air\n");
}

void METHANE19Kinetics::get_wdot_reduced(doublereal* wdot)
{
compute_reduced(thermo().pressure(), thermo().temperature(),thermo().massFractions(),wdot);
} 

void METHANE19Kinetics::compute_reduced(doublereal P, doublereal T, const doublereal* m_y,doublereal* wdot)
{
int i;
//pressure unit conversion for chemkin
P=P*10.0;
methane19_(&P,&T,&m_y[0],&wdot[0]);

//mol/mmol conversion cantera is in mmol
for (i=0;i<19;i++) {
	wdot[i]=1000.0*wdot[i];
}

}


}
