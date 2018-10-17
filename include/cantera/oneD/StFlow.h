/**
 * @file StFlow.h
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_STFLOW_H
#define CT_STFLOW_H

#include "cantera/transport/TransportBase.h"
#include "Domain1D.h"
#include "cantera/base/Array.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/numerics/funcs.h"

#include <cstdio>
#include <sstream>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace ctml;
using namespace std;

namespace Cantera
{
class MultiJac;

//------------------------------------------
//   constants
//------------------------------------------

// Offsets of solution components in the solution array.
const size_t c_offset_U = 0;    // axial velocity
const size_t c_offset_V = 1;    // strain rate
const size_t c_offset_T = 2;    // temperature
const size_t c_offset_L = 3;    // (1/r)dP/dr
const size_t c_offset_Y = 4;    // mass fractions

// Transport option flags
const int c_Mixav_Transport = 0;
const int c_Multi_Transport = 1;
const int c_Soret = 2;
const int c_AVBP_Transport = 3;

/**
 *  This class represents 1D flow domains that satisfy the one-dimensional
 *  similarity solution for chemically-reacting, axisymmetric, flows.
 *  @ingroup onedim
 */
class StFlow : public Domain1D
{
public:
    //--------------------------------
    // construction and destruction
    //--------------------------------

    //! Create a new flow domain.
    //! @param ph Object representing the gas phase. This object will be used
    //!     to evaluate all thermodynamic, kinetic, and transport properties.
    //! @param nsp Number of species.
    //! @param points Initial number of grid points
    StFlow(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! @name Problem Specification
    //! @{

    virtual void setupGrid(size_t n, const doublereal* z);

    thermo_t& phase() {
        return *m_thermo;
    }
    Kinetics& kinetics() {
        return *m_kin;
    }

    virtual void init() {
    }

    /**
     * Set the thermo manager. Note that the flow equations assume
     * the ideal gas equation.
     */
    void setThermo(IdealGasPhase& th) {
        m_thermo = &th;
    }

    //! Set the kinetics manager. The kinetics manager must
    void setKinetics(Kinetics& kin) {
        m_kin = &kin;
    }

    //! set the transport manager
    void setTransport(Transport& trans, bool withSoret = false);
    void enableSoret(bool withSoret);
    bool withSoret() const {
        return m_do_soret;
    }

    //! Set the pressure. Since the flow equations are for the limit of
    //! small Mach number, the pressure is very nearly constant
    //! throughout the flow.
    void setPressure(doublereal p) {
        m_press = p;
    }

    //! The current pressure [Pa].
    doublereal pressure() const {
        return m_press;
    }

    //! @deprecated unused
    virtual void setState(size_t point, const doublereal* state,
                          doublereal* x) {
        warn_deprecated("StFlow::setState");
        setTemperature(point, state[2]);
        for (size_t k = 0; k < m_nsp; k++) {
            setMassFraction(point, k, state[4+k]);
        }
    }

    //! Write the initial solution estimate into array x.
    virtual void _getInitialSoln(doublereal* x) {
        for (size_t j = 0; j < m_points; j++) {
            T(x,j) = m_thermo->temperature();
            m_thermo->getMassFractions(&Y(x, 0, j));
        }
    }

    virtual void _finalize(const doublereal* x);

    //! Sometimes it is desired to carry out the simulation using a specified
    //! temperature profile, rather than computing it by solving the energy
    //! equation. This method specifies this profile.
    void setFixedTempProfile(vector_fp& zfixed, vector_fp& tfixed) {
        m_zfix = zfixed;
        m_tfix = tfixed;
    }

    /*!
     * Set the temperature fixed point at grid point j, and disable the energy
     * equation so that the solution will be held to this value.
     */
    void setTemperature(size_t j, doublereal t) {
        m_fixedtemp[j] = t;
        m_do_energy[j] = false;
    }

    /*!
     * Set the mass fraction fixed point for species k at grid point j, and
     * disable the species equation so that the solution will be held to this
     * value. Note: in practice, the species are hardly ever held fixed.
     */
    void setMassFraction(size_t j, size_t k, doublereal y) {
        m_fixedy(k,j) = y;
        m_do_species[k] = true; // false;
    }

        /// Returns the mixture fraction at point j
        /// Added by P. Wolf - March 2010
        /// Thanks to Alireza Najafiyazdi
    doublereal getMixFrac(size_t j) {
                return m_zmixfrac[j];
        }

        /// Updates the mixture fraction at point j
        /// Added by P. Wolf - March 2010
        /// Thanks to Alireza Najafiyazdi
    void calcMixFrac(const doublereal* x,size_t j);
		


    //! The fixed temperature value at point j.
    doublereal T_fixed(size_t j) const {
        return m_fixedtemp[j];
    }

    //! The fixed mass fraction value of species k at point j.
    doublereal Y_fixed(size_t k, size_t j) const {
        return m_fixedy(k,j);
    }

    // @}

    virtual std::string componentName(size_t n) const;

    virtual size_t componentIndex(const std::string& name) const;

    //! Print the solution.
    virtual void showSolution(const doublereal* x);

    //! Save the current solution for this domain into an XML_Node
    /*!
     *  @param o    XML_Node to save the solution to.
     *  @param sol  Current value of the solution vector. The object will pick
     *              out which part of the solution vector pertains to this
     *              object.
     */
    virtual XML_Node& save(XML_Node& o, const doublereal* const sol);

    virtual void restore(const XML_Node& dom, doublereal* soln,
                         int loglevel);

    // overloaded in subclasses
    virtual std::string flowType() {
        return "<none>";
    }

    void solveEnergyEqn(size_t j=npos) {
        if (j == npos)
            for (size_t i = 0; i < m_points; i++) {
                m_do_energy[i] = true;
            }
        else {
            m_do_energy[j] = true;
        }
        m_refiner->setActive(0, true);
        m_refiner->setActive(1, true);
        m_refiner->setActive(2, true);
        needJacUpdate();
    }

    void fixTemperature(size_t j=npos) {
        if (j == npos)
            for (size_t i = 0; i < m_points; i++) {
                m_do_energy[i] = false;
            }
        else {
            m_do_energy[j] = false;
        }
        m_refiner->setActive(0, false);
        m_refiner->setActive(1, false);
        m_refiner->setActive(2, false);
        needJacUpdate();
    }

    bool doSpecies(size_t k) {
        return m_do_species[k];
    }
    bool doEnergy(size_t j) {
        return m_do_energy[j];
    }

    void solveSpecies(size_t k=npos) {
        if (k == npos) {
            for (size_t i = 0; i < m_nsp; i++) {
                m_do_species[i] = true;
            }
        } else {
            m_do_species[k] = true;
        }
        needJacUpdate();
    }

    void fixSpecies(size_t k=npos) {
        if (k == npos) {
            for (size_t i = 0; i < m_nsp; i++) {
                m_do_species[i] = false;
            }
        } else {
            m_do_species[k] = false;
        }
        needJacUpdate();
    }

    void integrateChem(doublereal* x,doublereal dt);

    //! Change the grid size. Called after grid refinement.
    void resize(size_t components, size_t points);

    virtual void setFixedPoint(int j0, doublereal t0) {}

    void setJac(MultiJac* jac);

    //! Set the gas object state to be consistent with the solution at point j.
    void setGas(const doublereal* x, size_t j);

    //! Set the gas state to be consistent with the solution at the midpoint
    //! between j and j + 1.
    void setGasAtMidpoint(const doublereal* x, size_t j);

    doublereal density(size_t j) const {
        return m_rho[j];
    }
        // P. Wolf
    void AVBPReadInputPea();

    //+++++++++++++++++++++++++++++++++++++++++++++++++++
    // B. Franzelli modifications to perform thickening
    //+++++++++++++++++++++++++++++++++++++++++++++++++++
    void AVBPAnalyticalOmegaZero();
    void AVBPReadInputChem();
    void AVBPcompute_local_thick(doublereal* x,size_t j);

    doublereal AVBPthick(size_t j) const {
            return avbp_thick[j];
    }
//   NOT implemented in this cantera version so must be carefull
    doublereal AVBPorder(size_t j, size_t i) const {
            return m_kin->reactantOrder(j,i);
    }

    virtual doublereal AVBPActivationTemp(doublereal *E);
//    doublereal AVBPActivationTemp(doublereal *E) const {
//        doublereal D[3];
//	if (m_kin->nReactions()==2) {
//	for (size_t n=0; n<m_kin->nReactions(); n++) {
//		this->m_kin->getArrheniusCoeffs(n, D);
//			// save activation temperature
//	            	E[n] = D[2];  
//		}
//            m_kin->getActivationEnergies(E);
//            return *E;
//    }
    /// getFwdRateOfProgress
    doublereal AVBPgetFwdRatesOfProgress(doublereal* x,size_t j, size_t avbp_reac_index) {
        setGas(x,j);
        doublereal *Fwdrate=new double[avbp_nreac];
        m_kin->getFwdRatesOfProgress(Fwdrate);
        doublereal temp = Fwdrate[avbp_reac_index];
        delete [] Fwdrate;
        return temp;
    }
    /// getNetRateOfProgress
    doublereal AVBPgetNetRatesOfProgress(doublereal* x,size_t j, size_t avbp_reac_index) {
        setGas(x,j);
        doublereal *Netrate=new double[avbp_nreac];
        m_kin->getNetRatesOfProgress(Netrate);
        doublereal temp = Netrate[avbp_reac_index];
        delete [] Netrate;
        return temp;
    }

    virtual bool fixed_mdot() {
        return true;
    }
    void setViscosityFlag(bool dovisc) {
        m_dovisc = dovisc;
    }

    /*!
     *  Evaluate the residual function for axisymmetric stagnation flow. If
     *  jpt is less than zero, the residual function is evaluated at all grid
     *  points. If jpt >= 0, then the residual function is only evaluated at
     *  grid points jpt-1, jpt, and jpt+1. This option is used to efficiently
     *  evaluate the Jacobian numerically.
     */
    virtual void eval(size_t j, doublereal* x, doublereal* r,
                      integer* mask, doublereal rdt);

    //! Evaluate all residual components at the right boundary.
    virtual void evalRightBoundary(doublereal* x, doublereal* res,
                                   integer* diag, doublereal rdt) = 0;

    //! Evaluate the residual corresponding to the continuity equation at all
    //! interior grid points.
    virtual void evalContinuity(size_t j, doublereal* x, doublereal* r,
                                integer* diag, doublereal rdt) = 0;

protected:
    doublereal component(const doublereal* x, size_t i, size_t j) const {
        return x[index(i,j)];
    }

    doublereal conc(const doublereal* x, size_t k,size_t j) const {
        return Y(x,k,j)*density(j)/m_wt[k];
    }

    doublereal cbar(const doublereal* x, size_t k, size_t j) const {
        return std::sqrt(8.0*GasConstant * T(x,j) / (Pi * m_wt[k]));
    }

    doublereal wdot(size_t k, size_t j) const {
        return m_wdot(k,j);
    }



/// write the net production rates at point j into array m_wdot
/// Modified by P. Wolf March 2010 - Thanks to Alireza N.
/// when PEA activated: production rates depends on mixture fraction at point j
/// WARNING: the following is valid only for fuel/AIR right now
    void getWdot(doublereal* x, size_t j) {
            		doublereal zloc;
            		doublereal phi;
          		doublereal C[3];

			setGas(x,j);		     // We'll need to set the gas before retrieving
			zloc =  getMixFrac(j);   // the Mixture Fraction and calculating wdot
			phi = zloc/(1. - zloc) * phi_cst;

			switch (avbp_ipea) {
				case 0:
					m_kin->getNetProductionRates(&m_wdot(0,j));
					break;
// ipea = 1, apply the PEA on fuel oxydation reaction only
				case 1:
					doublereal pea_f; //pea correction function
					doublereal A_sav; // Let's save the initial Pre-exponential
					// Now set the PEA correction function
					pea_f = 0.5*(1.+tanh((avbp_pea_coeffs[0]-phi)/avbp_pea_coeffs[1]))
					+ avbp_pea_coeffs[2]/4.*(1.+tanh((phi-avbp_pea_coeffs[3])/avbp_pea_coeffs[4]))
					* (1.+tanh((avbp_pea_coeffs[5]-phi)/avbp_pea_coeffs[6]));
					// Retrieve the Arrhenius coefficients
					this->m_kin->getArrheniusCoeffs(1, C);
	            			A_sav = C[0];	// Save A
					C[0] = C[0]*pea_f*pea_f;	// Apply PEA
	            	// Sets the modified Arrhenius and get Species net Production rates
					this->m_kin->setArrheniusCoeffs(1, C);
					m_kin->getNetProductionRates(&m_wdot(0,j));
            		// And restore the initial A
            				C[0] = A_sav;
					this->m_kin->setArrheniusCoeffs(1, C);
					break;
// ipea = 2, apply PEA on both reactions
				case 2:
//TEST 081014: differencie bien les cas
//                                        cout<<"Je suis dans le cas 2"<<endl;
					doublereal pea2_f[2]; //pea correction function
					doublereal A2_sav[2]; // Let's save the initial Pre-exponential
					// Now set the PEA correction function
					pea2_f[0] = 2./((1.+tanh((avbp_pea_coeffs[0]-phi)/avbp_pea_coeffs[1]))
						+ avbp_pea_coeffs[2]*(1.+tanh((phi-avbp_pea_coeffs[3])/avbp_pea_coeffs[4]))
						+ avbp_pea_coeffs[5]*(1.+tanh((phi-avbp_pea_coeffs[6])/avbp_pea_coeffs[7])));
				
					pea2_f[1] = 0.5*((1.+tanh((avbp_pea_coeffs[8]-phi)/avbp_pea_coeffs[9]))
					+ avbp_pea_coeffs[10]*(1.+tanh((phi-avbp_pea_coeffs[11])/avbp_pea_coeffs[12]))
					+ avbp_pea_coeffs[13]*(1.+tanh((phi-avbp_pea_coeffs[14])/avbp_pea_coeffs[15]))
					* (1.+tanh((avbp_pea_coeffs[16]-phi)/avbp_pea_coeffs[17])));
					
					// Retrieve the Arrhenius coefficients
					if (m_kin->nReactions()==2) {
						for (size_t n=0; n<m_kin->nReactions(); n++) {
							this->m_kin->getArrheniusCoeffs(n, C);
	            					A2_sav[n] = C[0];  // Save A
							C[0] = C[0]*pea2_f[n];
	            			// Sets the modified Arrhenius
							this->m_kin->setArrheniusCoeffs(n, C);
						}
	        			} else {
	        			throw CanteraError("getWdot","ipea = 2 but there ain't 2 reactions");
	          			}
	          			// get Species net Production rates
		            		m_kin->getNetProductionRates(&m_wdot(0,j));
					// And now, restore initial A - otherwise, we'll run into big troubles !!
					for (size_t n=0; n<m_kin->nReactions(); n++) {
						this->m_kin->getArrheniusCoeffs(n, C);
	           				C[0] = A2_sav[n];
						this->m_kin->setArrheniusCoeffs(n, C);
	        			}
					break;
// ipea = 3, apply PEA on both reactions (for CM2 only)
				case 3:
					doublereal pea3_f[2]; //pea correction function
					doublereal A3_sav[2]; // Let's save the initial Pre-exponential
					// Now set the PEA correction function
					if (phi>1.0) { 
                                          //pea3_f[0] = max(avbp_pea_coeffs[0]*phi*phi + avbp_pea_coeffs[1]*phi + avbp_pea_coeffs[2],0.0); 
                                          pea3_f[0] = avbp_pea_coeffs[0]*phi*phi + avbp_pea_coeffs[1]*phi + avbp_pea_coeffs[2]; 
					} else { 
                                          pea3_f[0] = 1; 
                                        }
					pea3_f[1] = pea3_f[0];
                                        //cout<<endl<<"check coeff "<<avbp_pea_coeffs[0]<<avbp_pea_coeffs[1]<<avbp_pea_coeffs[2]<<endl;
					//cout<<endl<<"check PEA 3 "<<phi<<" =>  "<<pea3_f[0];
					// Retrieve the Arrhenius coefficients
					if (m_kin->nReactions()==2) {
						for (size_t n=0; n<m_kin->nReactions(); n++) {
							this->m_kin->getArrheniusCoeffs(n, C);
	            					A3_sav[n] = C[0];  // Save A
							C[0] = C[0]*pea3_f[n];
	            			// Sets the modified Arrhenius
							this->m_kin->setArrheniusCoeffs(n, C);
						}
	        			} else {
	        				throw CanteraError("getWdot","ipea = 2 but there ain't 2 reactions");
	          			}
	          			// get Species net Production rates
		           		 m_kin->getNetProductionRates(&m_wdot(0,j));
					// And now, restore initial A - otherwise, we'll run into big troubles !!
					for (size_t n=0; n<m_kin->nReactions(); n++) {
						this->m_kin->getArrheniusCoeffs(n, C);
	           				C[0] = A3_sav[n];
						this->m_kin->setArrheniusCoeffs(n, C);
	        			}
					break;

				default:
					throw CanteraError("getWdot","ipea is set to a weird value");
					break;
			}			
    }

    /**
     * Update the thermodynamic properties from point j0 to point j1
     * (inclusive), based on solution x.
     */
    void updateThermo(const doublereal* x, size_t j0, size_t j1) {
        for (size_t j = j0; j <= j1; j++) {
            setGas(x,j);
            m_rho[j] = m_thermo->density();
            m_wtm[j] = m_thermo->meanMolecularWeight();
            m_cp[j]  = m_thermo->cp_mass();
        }
    }

    //--------------------------------
    // central-differenced derivatives
    //--------------------------------

    doublereal cdif2(const doublereal* x, size_t n, size_t j,
                     const doublereal* f) const {
        doublereal c1 = (f[j] + f[j-1])*(x[index(n,j)] - x[index(n,j-1)]);
        doublereal c2 = (f[j+1] + f[j])*(x[index(n,j+1)] - x[index(n,j)]);
        return (c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }


    //! @name Solution components
    //! @{

    doublereal T(const doublereal* x, size_t j) const {
        return x[index(c_offset_T, j)];
    }
    doublereal& T(doublereal* x, size_t j) {
        return x[index(c_offset_T, j)];
    }
    doublereal T_prev(size_t j) const {
        return prevSoln(c_offset_T, j);
    }

    doublereal rho_u(const doublereal* x, size_t j) const {
        return m_rho[j]*x[index(c_offset_U, j)];
    }

    doublereal u(const doublereal* x, size_t j) const {
        return x[index(c_offset_U, j)];
    }

    doublereal V(const doublereal* x, size_t j) const {
        return x[index(c_offset_V, j)];
    }
    doublereal V_prev(size_t j) const {
        return prevSoln(c_offset_V, j);
    }

    doublereal lambda(const doublereal* x, size_t j) const {
        return x[index(c_offset_L, j)];
    }

    doublereal Y(const doublereal* x, size_t k, size_t j) const {
        return x[index(c_offset_Y + k, j)];
    }

    doublereal& Y(doublereal* x, size_t k, size_t j) {
        return x[index(c_offset_Y + k, j)];
    }

    doublereal Y_prev(size_t k, size_t j) const {
        return prevSoln(c_offset_Y + k, j);
    }

    doublereal X(const doublereal* x, size_t k, size_t j) const {
        return m_wtm[j]*Y(x,k,j)/m_wt[k];
    }

    doublereal flux(size_t k, size_t j) const {
        return m_flux(k, j);
    }
    //! @}

    //! @name convective spatial derivatives.
    //! These use upwind differencing, assuming u(z) is negative
    //! @{
    doublereal dVdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (V(x,jloc) - V(x,jloc-1))/m_dz[jloc-1];
    }

    doublereal dYdz(const doublereal* x, size_t k, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (Y(x,k,jloc) - Y(x,k,jloc-1))/m_dz[jloc-1];
    }

    doublereal dTdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (T(x,jloc) - T(x,jloc-1))/m_dz[jloc-1];
    }
    //! @}

    doublereal shear(const doublereal* x, size_t j) const {
        doublereal c1 = m_visc[j-1]*(V(x,j) - V(x,j-1));
        doublereal c2 = m_visc[j]*(V(x,j+1) - V(x,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    doublereal divHeatFlux(const doublereal* x, size_t j) const {
        doublereal c1 = m_tcon[j-1]*(T(x,j) - T(x,j-1));
        doublereal c2 = m_tcon[j]*(T(x,j+1) - T(x,j));
        return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    doublereal AVBPdivHeatFlux(const doublereal* x, size_t j) const {
            doublereal c1 = avbp_thick[j-1]*m_tcon[j-1]*(T(x,j) - T(x,j-1));
            doublereal c2 = avbp_thick[j]*m_tcon[j]*(T(x,j+1) - T(x,j));
            return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    size_t mindex(size_t k, size_t j, size_t m) {
        return m*m_nsp*m_nsp + m_nsp*j + k;
    }

    //! Update the diffusive mass fluxes.
    void updateDiffFluxes(const doublereal* x, size_t j0, size_t j1);

    //---------------------------------------------------------
    //             member data
    //---------------------------------------------------------

    // inlet
    doublereal m_inlet_u;
    doublereal m_inlet_V;
    doublereal m_inlet_T;
    doublereal m_rho_inlet;
    vector_fp m_yin;

    // surface
    doublereal m_surface_T;

    doublereal m_press;        // pressure

    // grid parameters
    vector_fp m_dz;
    //vector_fp m_z;

    // mixture thermo properties
    vector_fp m_rho;
    vector_fp m_wtm;
// P. Wolf
    vector_fp m_zmixfrac;

    // species thermo properties
    vector_fp m_wt;
    vector_fp m_cp;
    vector_fp m_enth;

        // Pierre WOLF: added parameters for input_pea.dat
    size_t avbp_ipea;
    std::string avbp_fuel;
    vector_fp avbp_pea_coeffs;
    doublereal phi_cst;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++
    // B. Franzelli modifications to perform thickening
    // input_chem parameters
    //+++++++++++++++++++++++++++++++++++++++++++++++++++
    size_t avbp_ithick;
    doublereal avbp_fthick;
    vector_fp avbp_thick;
    size_t avbp_iclip;
    doublereal avbp_delta_chem;
    doublereal avbp_betaprime;
    doublereal avbp_ta;
    doublereal avbp_tin;
    doublereal avbp_tout;
    doublereal avbp_sratio;
    doublereal avbp_omega_zero;
    doublereal avbp_omega;
    doublereal theta_F;
    size_t avbp_iomega;
    size_t avbp_reac_index;
    size_t avbp_fuel_index;
    size_t avbp_oxy_index;
    size_t avbp_nreac;




    // transport properties
    vector_fp m_visc;
    vector_fp m_tcon;
    vector_fp m_diff;
    vector_fp m_multidiff;
    Array2D m_dthermal;
    Array2D m_flux;

    // production rates
    Array2D m_wdot;
    vector_fp m_surfdot;

    size_t m_nsp;

    IdealGasPhase* m_thermo;
    Kinetics* m_kin;
    Transport* m_trans;

    MultiJac* m_jac;

    bool m_ok;

    // flags
    std::vector<bool> m_do_energy;
    bool m_do_soret;
    std::vector<bool> m_do_species;
    int m_transport_option;

    // solution estimate
    //vector_fp m_zest;
    //Array2D   m_yest;

    // fixed T and Y values
    Array2D   m_fixedy;
    vector_fp m_fixedtemp;
    vector_fp m_zfix;
    vector_fp m_tfix;

    bool m_dovisc;
    doublereal m_efctr;

    //! Update the transport properties at grid points in the range from `j0`
    //! to `j1`, based on solution `x`.
    void updateTransport(doublereal* x, size_t j0, size_t j1);

private:
    vector_fp m_ybar;
};

/**
 * A class for axisymmetric stagnation flows.
 * @ingroup onedim
 */
class AxiStagnFlow : public StFlow
{
public:
    AxiStagnFlow(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1) :
        StFlow(ph, nsp, points) {
        m_dovisc = true;
    }

    virtual void evalRightBoundary(doublereal* x, doublereal* res,
                                   integer* diag, doublereal rdt);
    virtual void evalContinuity(size_t j, doublereal* x, doublereal* r,
                                integer* diag, doublereal rdt);

    virtual std::string flowType() {
        return "Axisymmetric Stagnation";
    }
};

/**
 * A class for freely-propagating premixed flames.
 * @ingroup onedim
 */
class FreeFlame : public StFlow
{
public:
    FreeFlame(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1) :
        StFlow(ph, nsp, points) {
        m_dovisc = false;
        setID("flame");
    }

    virtual void evalRightBoundary(doublereal* x, doublereal* res,
                                   integer* diag, doublereal rdt);
    virtual void evalContinuity(size_t j, doublereal* x, doublereal* r,
                                integer* diag, doublereal rdt);

    virtual std::string flowType() {
        return "Free Flame";
    }
    virtual bool fixed_mdot() {
        return false;
    }
};

/**
 * Import a previous solution to use as an initial estimate. The
 * previous solution may have been computed using a different
 * reaction mechanism. Species in the old and new mechanisms are
 * matched by name, and any species in the new mechanism that were
 * not in the old one are set to zero. The new solution is created
 * with the same number of grid points as in the old solution.
 * @deprecated Conversion is handled automatically when loading from XML
 */
void importSolution(size_t points, doublereal* oldSoln, IdealGasPhase& oldmech,
                    size_t size_new, doublereal* newSoln, IdealGasPhase& newmech);

}

#endif
