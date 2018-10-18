/**
 * @file ETHANOL28Kinetics.h
 *
 * @ingroup chemkinetics
 */

// Author T. Jaravel

#ifndef CT_ETHANOL28_KINETICS_H
#define CT_ETHANOL28_KINETICS_H

#include "GasKinetics.h"

namespace Cantera
{
const int cETHANOL28Kinetics = 155;

/**
 *  Kinetics manager implementing reaction mechanism GRI-Mech 3.0
 *  @deprecated
 */
class ETHANOL28Kinetics : public GasKinetics
{
public:
    /// Default constructor.
    ETHANOL28Kinetics(thermo_t* th=0);

    virtual int type() const {
        return cETHANOL28Kinetics;
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
