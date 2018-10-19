/**
 *
 * @file METHANE13Kinetics.h
 *
 * @ingroup chemkinetics
 */
// Author T. Jaravel


#ifndef CT_METHANE13_KINETICS_H
#define CT_METHANE13_KINETICS_H

#include "GasKinetics.h"

namespace Cantera {

    const int cMETHANE13Kinetics = 150;

    /**
     *  Kinetics manager implementing reaction mechanism METHANE
     */    
    class METHANE13Kinetics : public GasKinetics {

    public:

        /// Default constructor.
        METHANE13Kinetics(thermo_t* th=0);

        virtual int type() const { return cMETHANE13Kinetics; }

        virtual void getNetProductionRates(doublereal* net) {
            get_wdot_reduced(net);
        }
      
    private:
        void get_wdot_reduced(doublereal* wdot);
        void compute_reduced(doublereal P, doublereal T,const doublereal* m_y,doublereal* wdot);
    };
}

#endif
