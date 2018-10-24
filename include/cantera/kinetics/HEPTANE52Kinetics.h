/**
 * @file HEPTANE52Kinetics.h
 *
 * @ingroup chemkinetics
 */

// Author T. Jaravel

#ifndef CT_HEPTANE52_KINETICS_H
#define CT_HEPTANE52_KINETICS_H

#include "GasKinetics.h"

namespace Cantera
{
const int cHEPTANE52Kinetics = 105;

/**
 *  Kinetics manager implementing reaction mechanism GRI-Mech 3.0
 *  @deprecated
 */
class HEPTANE52Kinetics : public GasKinetics
{
public:
    /// Default constructor.
    HEPTANE52Kinetics(thermo_t* th=0);

    virtual int type() const {
        return cHEPTANE52Kinetics;
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
