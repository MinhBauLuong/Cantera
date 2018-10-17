/**
 *  @file DME30Kinetics.cpp 
 *
 * @ingroup chemkinetics
 */

// Author: S. Desai

#include "cantera/kinetics/DME30Kinetics.h"

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Enhanced3BConc.h"
#include "cantera/kinetics/ThirdBodyMgr.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/thermo/IdealGasPhase.h"
#include <iostream>

using namespace std;


namespace Cantera {
  //Fortran External Routine
  extern "C" {
    void dme30_(doublereal* P, doublereal* T, const doublereal* m_y,doublereal* wdot);
          }
    /**
     * Construct an empty reaction mechanism.
     */
      
    DME30Kinetics::
    DME30Kinetics(thermo_t* th) : GasKinetics(th) {
    //printf("Warning: You are using a 13 species reduced scheme for methane air oxidation");
    }


   void DME30Kinetics::get_wdot_reduced(double* wdot){
        compute_reduced(thermo().pressure(), thermo().temperature(),thermo().massFractions(),wdot);
    } 


   void DME30Kinetics::compute_reduced(doublereal P, doublereal T, const doublereal* m_y,doublereal* wdot){
      int i;

      //pressure unit conversion for chemkin
      P=P*10.0;
      dme30_(&P,&T,&m_y[0],&wdot[0]);

      //mol/mmol conversion cantera is in mmol
      for (i=0;i<30;i++) {
              wdot[i]=1000.0*wdot[i];
       }

   }



}








