/**
 * @file TESTKinetics.h
 *
 * @ingroup chemkinetics
 */

// Author T. Jaravel

#ifndef CT_TEST_KINETICS_H
#define CT_TEST_KINETICS_H

#include "GasKinetics.h"

namespace Cantera
{
const int cTESTKinetics = 155;

/**
 *  Kinetics manager implementing reaction mechanism GRI-Mech 3.0
 *  @deprecated
 */
class TESTKinetics : public GasKinetics
{
public:
    /// Default constructor.
    TESTKinetics(thermo_t* th=0);

    virtual int type() const {
        return cTESTKinetics;
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
