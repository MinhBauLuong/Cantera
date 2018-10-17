/**
 *
 * @file DME30Kinetics.h
 *
 * @ingroup chemkinetics
 */
// Author S. Desai


#ifndef CT_DME30_KINETICS_H
#define CT_DME30_KINETICS_H

#include "GasKinetics.h"

namespace Cantera {

    const int cDME30Kinetics = 111;

    /**
     *  Kinetics manager implementing reaction mechanism DME30
     */    
    class DME30Kinetics : public GasKinetics {

    public:

        /// Default constructor.
        DME30Kinetics(thermo_t* th=0);

        virtual int type() const { return cDME30Kinetics; }

        virtual void getNetProductionRates(doublereal* net) {
            get_wdot_reduced(net);
        }
      
    private:
        void get_wdot_reduced(doublereal* wdot);
        void compute_reduced(doublereal P, doublereal T,const doublereal* m_y,doublereal* wdot);
    };
}

#endif
