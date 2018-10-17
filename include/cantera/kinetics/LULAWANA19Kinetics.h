/**
 * @file LULAWANA19Kinetics.h
 *
 * @ingroup chemkinetics
 */

// Author T. Jaravel

#ifndef CT_LULAWANA19_KINETICS_H
#define CT_LULAWANA19_KINETICS_H

#include "GasKinetics.h"

namespace Cantera
{
const int cLULAWANA19Kinetics = 101;

/**
 *  Kinetics manager implementing reaction mechanism GRI-Mech 3.0
 *  @deprecated
 */
class LULAWANA19Kinetics : public GasKinetics
{
public:
    /// Default constructor.
    LULAWANA19Kinetics(thermo_t* th=0);

    virtual int type() const {
        return cLULAWANA19Kinetics;
    }

    virtual void getNetProductionRates(doublereal* net) {
          get_wdot_reduced(net);
    }

private:
    void get_wdot_reduced(doublereal* wdot);
    void compute_reduced(doublereal P, doublereal T,const doublereal* m_y,doublereal* wdot);
};
}


#endif
