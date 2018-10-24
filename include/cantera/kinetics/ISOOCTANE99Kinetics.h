/**
 * @file ISOOCTANE99Kinetics.h
 *
 * @ingroup chemkinetics
 */

// Author T. Jaravel

#ifndef CT_ISOOCTANE99_KINETICS_H
#define CT_ISOOCTANE99_KINETICS_H

#include "GasKinetics.h"

namespace Cantera
{
const int cISOOCTANE99Kinetics = 109;

/**
 *  Kinetics manager implementing reaction mechanism GRI-Mech 3.0
 *  @deprecated
 */
class ISOOCTANE99Kinetics : public GasKinetics
{
public:
    /// Default constructor.
    ISOOCTANE99Kinetics(thermo_t* th=0);

    virtual int type() const {
        return cISOOCTANE99Kinetics;
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
