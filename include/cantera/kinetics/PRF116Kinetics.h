/**
 * @file PRF116Kinetics.h
 *
 * @ingroup chemkinetics
 */

// Author T. Jaravel

#ifndef CT_PRF116_KINETICS_H
#define CT_PRF116_KINETICS_H

#include "GasKinetics.h"

namespace Cantera
{
const int cPRF116Kinetics = 158;

/**
 *  Kinetics manager implementing reaction mechanism GRI-Mech 3.0
 *  @deprecated
 */
class PRF116Kinetics : public GasKinetics
{
public:
    /// Default constructor.
    PRF116Kinetics(thermo_t* th=0);

    virtual int type() const {
        return cPRF116Kinetics;
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
