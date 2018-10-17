/**
 *  @file AVBPTransport.cpp
 *  Simplified AVBP transport properties for ideal gas mixtures.
 */

/* $Author: franzelli $
 * $Revision: felden $
 * $Date: 7/10/2014 $
 */
// copyright 2001 California Institute of Technology

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/AVBPTransport.h"
#include "cantera/base/utilities.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/IdealGasPhase.h"

#include <sstream>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

namespace Cantera
{

AVBPTransport::AVBPTransport() :
    m_nsp(0),
    m_temp(-1.0),
    m_logt(0.0)
{
}

bool AVBPTransport::initGas(GasTransportParams& tr)
{
    // constant substance attributes
    m_thermo = tr.thermo;
    m_nsp   = m_thermo->nSpecies();

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(),
         m_thermo->molecularWeights().end(), m_mw.begin());

    // copy polynomials and parameters into local storage
    m_poly       = tr.poly;
    m_visccoeffs = tr.visccoeffs;
    m_condcoeffs = tr.condcoeffs;
    m_diffcoeffs = tr.diffcoeffs;

    m_zrot       = tr.zrot;
    m_crot       = tr.crot;
    m_epsilon    = tr.epsilon;
    m_mode       = tr.mode_;
    m_diam       = tr.diam;
    m_eps        = tr.eps;
    m_alpha      = tr.alpha;
    m_dipoleDiag.resize(m_nsp);
    for (int i = 0; i < m_nsp; i++) {
        m_dipoleDiag[i] = tr.dipole(i,i);
    }

    m_phi.resize(m_nsp, m_nsp, 0.0);
    m_wratjk.resize(m_nsp, m_nsp, 0.0);
    m_wratkj1.resize(m_nsp, m_nsp, 0.0);
    int j, k;
    for (j = 0; j < m_nsp; j++)
        for (k = j; k < m_nsp; k++) {
            m_wratjk(j,k) = sqrt(m_mw[j]/m_mw[k]);
            m_wratjk(k,j) = sqrt(m_wratjk(j,k));
            m_wratkj1(j,k) = sqrt(1.0 + m_mw[k]/m_mw[j]);
        }

    m_polytempvec.resize(5);
    m_visc.resize(m_nsp);
    m_sqvisc.resize(m_nsp);
    m_cond.resize(m_nsp);
    m_bdiff.resize(m_nsp, m_nsp);

    m_molefracs.resize(m_nsp);
    m_spwork.resize(m_nsp);

    // set flags all false
    m_viscmix_ok = false;
    m_viscwt_ok = false;
    m_spvisc_ok = false;
    m_spcond_ok = false;
    m_condmix_ok = false;
    m_spcond_ok = false;
    m_diffmix_ok = false;
    m_abc_ok = false;

//AVBP 071014 INPUT THERMO
    //lecture du fichier input_thermo et stockage des donnees en char
    ifstream avbp_thermo;
    char dummy[200];
//
    char avbp_Prandtll[13];
    double avbp_Prandtl_d;
//
    char avbp_mu00[13];
    double avbp_mu0_d;
//
    char avbp_T00[13];
    double avbp_T0_d;
//
    char avbp_betaa[13];
    double avbp_beta_d;

    avbp_thermo.open("input_thermo.dat");                                                              
    if (!avbp_thermo){                                                                                 
          cout<<"MANNAGGIA: can not find the 'input_thermo.dat' file!!!";                                  
          return -1;
        }                                                                            
    avbp_thermo>>avbp_mu00;                                                                                                                                            
    avbp_thermo.getline(dummy,200,'\n');
    avbp_thermo>>avbp_T00;                                                                                                                                            
    avbp_thermo.getline(dummy,200,'\n');
    avbp_thermo>>avbp_betaa;                                                                                                                                            
    avbp_thermo.getline(dummy,200,'\n');
    avbp_thermo>>avbp_Prandtll;                                                                                                                                            
    avbp_thermo.getline(dummy,200,'\n');
    avbp_thermo.close();

    //conversion donnees en float lisible cpp tjs avec des char
    for (int i=0; i<12; i++){
        if (avbp_mu00[i]=='d'){
                avbp_mu00[i]='e';
		cout<<"WARNING fichier input_thermo illisible ... mais je gère"<<endl;
        }
    }
    for (int i=0; i<12; i++){
        if (avbp_T00[i]=='d'){
                avbp_T00[i]='e';
		cout<<"WARNING fichier input_thermo illisible ... mais je gère"<<endl;
        }
    }
    for (int i=0; i<12; i++){
        if (avbp_betaa[i]=='d'){
                avbp_betaa[i]='e';
		cout<<"WARNING fichier input_thermo illisible ... mais je gère"<<endl;
        }
    }
    for (int i=0; i<12; i++){
        if (avbp_Prandtll[i]=='d'){
                avbp_Prandtll[i]='e';
		cout<<"WARNING fichier input_thermo illisible ... mais je gère"<<endl;
        }
    }

	//conversion char en double
    avbp_mu0_d=atof(avbp_mu00);
    avbp_T0_d=atof(avbp_T00);
    avbp_beta_d=atof(avbp_betaa);
    avbp_Prandtl_d=atof(avbp_Prandtll);

	//stockage bonnes donnees
    avbp_mu0=avbp_mu0_d;
    avbp_T0=avbp_T0_d;
    avbp_beta=avbp_beta_d;
    avbp_Prandtl=avbp_Prandtl_d;

//AVBP 071014 INPUT PREMIX
    //lecture du fichier input_premix et stockage des donnees en char
    ifstream avbp_premix;
    char avbp_Schh[13];
    double avbp_Sch_d;
    avbp_premix.open("input_premix.dat");
    if (!avbp_premix){
      cout<<"MANNAGGIA: can not find the 'input_premix.dat' file!!";
      return -1;
    }
    //ignorer les 4 prem lignes
    avbp_Le.resize(m_nsp);
    avbp_Sch.resize(m_nsp);
    for(int k=0;k<4;k++){  
      avbp_premix.ignore(80,'\n');
    }  
    //lecture du fichier input_premix et modif avant stockage
    for(int i=0;i<m_nsp;i++){
      avbp_premix>>dummy;
      avbp_premix>>avbp_Schh;
      for (int i=0; i<12; i++){
	      if (avbp_Schh[i]=='d'){
                  avbp_Schh[i]='e';
		  cout<<"WARNING fichier input_premix illisible ... mais je gère"<<endl;
	      }
      }
      avbp_Sch_d=atof(avbp_Schh);
      avbp_Sch[i]=avbp_Sch_d;
      avbp_Le[i]= avbp_Sch[i] / avbp_Prandtl;
      avbp_premix.ignore(80,'\n'); 
    }
    avbp_premix.close();

//AVBP 071014 INITIALISATION
    cout<<endl<<"++++++++++AVBP Transports is used:++++++++++++"<<endl;
    cout<<"-Reference viscosity: "<<avbp_mu0<<endl;
    cout<<"-Reference temperature: "<<avbp_T0<<endl;
    if(avbp_beta<0){
      cout<<"-Exponent for power law: "<<-avbp_beta<<endl;
    }
    else if(avbp_beta>0){
       cout<<"-Exponent for Second Sutherland Constant: "<<avbp_beta<<endl;
    }
    else{
       cout<<"The molecular viscosity used is temperature indipendent and equal to: "<<avbp_mu0;
    }
    cout<<"-Prandtl number: "<<avbp_Prandtl<<endl;
    cout<<"-Schmidt and Lewis number for: "<<endl;
    for(int i=0;i<m_nsp;i++){
      cout<<"*Specie "<<i<<": "<<avbp_Sch[i]<<" "<<avbp_Le[i]<<endl;
    }

    doublereal m_cp = 0.0;
    doublereal m_density = 0.0;
    m_cp = m_thermo->cp_mass();
    m_density = m_thermo->density();
    if(avbp_beta<0){
      doublereal avbp_absbeta;
      avbp_absbeta = - avbp_beta;            
      m_viscmix = avbp_mu0 * pow(m_temp/avbp_T0,avbp_absbeta);
    }
    else if(avbp_beta>0){
      doublereal coeff;
      coeff= (avbp_T0 + avbp_beta) / pow(avbp_T0,1.5); 
      m_viscmix = avbp_mu0 * coeff * pow(m_temp, 1.5) / (m_temp + avbp_beta);  
    }
    else{
    m_viscmix = avbp_mu0;
    }
    m_lambda = m_viscmix * m_cp / avbp_Prandtl;


    return true;
}


doublereal AVBPTransport::viscosity()
{
    update_T();
    update_C();

    if (m_viscmix_ok) {
        return m_viscmix;
    }

    doublereal vismix = 0.0;
    // AVBP Transport Properties
    if(avbp_beta<0){
          doublereal avbp_absbeta;
          avbp_absbeta = - avbp_beta;
          vismix = avbp_mu0 * pow(m_temp/avbp_T0,avbp_absbeta);
    }
    else if(avbp_beta>0){
          doublereal coeff;
          coeff= (avbp_T0 + avbp_beta) / pow(avbp_T0,1.5);
          vismix = avbp_mu0 * coeff * pow(m_temp, 1.5) / (m_temp + avbp_beta);
    }
    else{
        vismix = avbp_mu0;
    }   
     
    m_viscmix = vismix;
    return m_viscmix;
}

void AVBPTransport::getBinaryDiffCoeffs(const size_t ld, doublereal* const d)
{
    int i,j;

    update_T();

    // if necessary, evaluate the binary diffusion coefficents
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    doublereal rp = 1.0/pressure_ig();
    for (i = 0; i < m_nsp; i++)
        for (j = 0; j < m_nsp; j++) {
            d[ld*j + i] = rp * m_bdiff(i,j);
        }
}

void AVBPTransport::getMobilities(doublereal* const mobil)
{
    int k;
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k] * m_thermo->charge(k);
    }
}

doublereal AVBPTransport::thermalConductivity()
{
    int k;
    doublereal lambda = 0.0;

    update_T();
    update_C();

    // update m_cond and m_phi if necessary
    if (!m_spcond_ok) {
        updateCond_T();
    }
        // AVBP transport properties 
    doublereal m_cp = 0.0;
    doublereal m_vismix = 0.0;
    m_cp = m_thermo->cp_mass();
    
    if(avbp_beta<0){
          doublereal avbp_absbeta;
          avbp_absbeta = - avbp_beta; 
          m_vismix = avbp_mu0 * pow(m_temp/avbp_T0,avbp_absbeta);
    }else if(avbp_beta>0){   
          doublereal coeff;
          coeff = (avbp_T0 + avbp_beta) / pow(avbp_T0,1.5); 
          m_vismix = avbp_mu0 * coeff * pow (m_temp, 1.5) / (m_temp + avbp_beta);
    }else{
          m_vismix = avbp_mu0;
    }

    m_lambda = m_vismix * m_cp / avbp_Prandtl;
    return m_lambda;
}

void AVBPTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    int k;
    for (k = 0; k < m_nsp; k++) {
        dt[k] = 0.0;
    }
}

void AVBPTransport::getSpeciesFluxes(size_t ndim,
                                      const doublereal* const grad_T,
                                      size_t ldx, const doublereal* const grad_X,
                                      size_t ldf, doublereal* const fluxes)
{
    size_t n = 0;
    int k;

    update_T();
    update_C();

    getMixDiffCoeffs(DATA_PTR(m_spwork));

    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* y  = m_thermo->massFractions();
    doublereal rhon = m_thermo->molarDensity();

    vector_fp sum(ndim,0.0);

    doublereal correction=0.0;
    // grab 2nd (summation) term -- still need to multiply by mass fraction (\rho_s / \rho)
    for (k = 0; k < m_nsp; k++) {
        correction += rhon * mw[k] * m_spwork[k] * grad_X[n*ldx + k];
    }

    for (n = 0; n < ndim; n++) {
        for (k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * grad_X[n*ldx + k] + y[k]*correction;
            sum[n] += fluxes[n*ldf + k];
        }
    }
    // add correction flux to enforce sum to zero
    for (n = 0; n < ndim; n++) {
        for (k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] -= y[k]*sum[n];
        }
    }
}

void AVBPTransport::getMixDiffCoeffs(doublereal* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    int k, j;
    doublereal mmw = m_thermo->meanMolecularWeight();
    doublereal sumxw = 0.0, sum2;
    doublereal p = pressure_ig();
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (k = 0; k < m_nsp; k++) {
            sumxw += m_molefracs[k] * m_mw[k];
        }
        for (k = 0; k < m_nsp; k++) {
            sum2 = 0.0;
            for (j = 0; j < m_nsp; j++) {
                if (j != k) {
                    sum2 += m_molefracs[j] / m_bdiff(j,k);
                }
            }
            if (sum2 <= 0.0) {
                d[k] = m_bdiff(k,k) / p;
            } else {
                d[k] = (sumxw - m_molefracs[k] * m_mw[k])/(p * mmw * sum2);
            }
        }
    }
         // AVBP Transport Properties 
    doublereal m_cp = 0.0;
    doublereal m_density = 0.0;
    m_cp = m_thermo->cp_mass();
    m_density = m_thermo->density();
 
    doublereal m_vismix = 0.0;
    if(avbp_beta<0){
          doublereal avbp_absbeta;
          avbp_absbeta = - avbp_beta;
          m_vismix = avbp_mu0 * pow(m_temp/avbp_T0,avbp_absbeta);
    }else if(avbp_beta>0){      
          doublereal coeff;
          coeff = (avbp_T0 + avbp_beta) / pow(avbp_T0,1.5);
          m_vismix = avbp_mu0 * coeff * pow(m_temp,1.5)/ (m_temp + avbp_beta);
    }else{
          m_vismix = avbp_mu0;
    }
 
    for (int k=0; k<m_nsp; k++){
          d[k]=m_vismix / m_density / avbp_Sch[k];
    }
}

/**
 *  @internal This is called whenever a transport property is
 *  requested from ThermoSubstance if the temperature has changed
 *  since the last call to update_T.
 */
void AVBPTransport::update_T()
{
    doublereal t = m_thermo->temperature();
    if (t == m_temp) {
        return;
    }
    if (t <= 0.0) {
        throw CanteraError("AVBPTransport::update_T",
                           "negative temperature "+fp2str(t));
    }
    m_temp = t;
    m_logt = log(m_temp);
    m_kbt = Boltzmann * m_temp;
    m_sqrt_t = sqrt(m_temp);
    m_t14 = sqrt(m_sqrt_t);
    m_t32 = m_temp * m_sqrt_t;
    m_sqrt_kbt = sqrt(Boltzmann*m_temp);

    // compute powers of log(T)
    m_polytempvec[0] = 1.0;
    m_polytempvec[1] = m_logt;
    m_polytempvec[2] = m_logt*m_logt;
    m_polytempvec[3] = m_logt*m_logt*m_logt;
    m_polytempvec[4] = m_logt*m_logt*m_logt*m_logt;

    // temperature has changed, so polynomial fits will need to be redone.
    m_viscmix_ok = false;
    m_spvisc_ok = false;
    m_viscwt_ok = false;
    m_spcond_ok = false;
    m_diffmix_ok = false;
    m_bindiff_ok = false;
    m_abc_ok  = false;
    m_condmix_ok = false;
}

void AVBPTransport::update_C()
{
    // signal that concentration-dependent quantities will need to
    // be recomputed before use, and update the local mole
    // fractions.

    m_viscmix_ok = false;
    m_diffmix_ok = false;
    m_condmix_ok = false;

    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

    // add an offset to avoid a pure species condition
    int k;
    for (k = 0; k < m_nsp; k++) {
        m_molefracs[k] = std::max(Tiny, m_molefracs[k]);
    }
}

void AVBPTransport::updateCond_T()
{
    if (m_mode == CK_Mode) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = exp(dot4(m_polytempvec, m_condcoeffs[k]));
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = m_sqrt_t * dot5(m_polytempvec, m_condcoeffs[k]);
        }
    }
    m_spcond_ok = true;
    m_condmix_ok = false;
}

void AVBPTransport::updateDiff_T()
{
    // evaluate binary diffusion coefficients at unit pressure
    int i,j;
    int ic = 0;
    if (m_mode == CK_Mode) {
        for (i = 0; i < m_nsp; i++) {
            for (j = i; j < m_nsp; j++) {
                m_bdiff(i,j) = exp(dot4(m_polytempvec, m_diffcoeffs[ic]));
                m_bdiff(j,i) = m_bdiff(i,j);
                ic++;
            }
        }
    } else {
        for (i = 0; i < m_nsp; i++) {
            for (j = i; j < m_nsp; j++) {
                m_bdiff(i,j) = m_temp * m_sqrt_t*dot5(m_polytempvec,
                                                      m_diffcoeffs[ic]);
                m_bdiff(j,i) = m_bdiff(i,j);
                ic++;
            }
        }
    }

    m_bindiff_ok = true;
    m_diffmix_ok = false;
}

void AVBPTransport::updateSpeciesViscosities()
{
    update_T();
    if (m_mode == CK_Mode) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_visc[k] = exp(dot4(m_polytempvec, m_visccoeffs[k]));
            m_sqvisc[k] = sqrt(m_visc[k]);
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            // the polynomial fit is done for sqrt(visc/sqrt(T))
            m_sqvisc[k] = m_t14 * dot5(m_polytempvec, m_visccoeffs[k]);
            m_visc[k] = (m_sqvisc[k] * m_sqvisc[k]);
        }
    }
    m_spvisc_ok = true;
}

void AVBPTransport::updateViscosity_T()
{
    doublereal vratiokj, wratiojk, factor1;

    if (!m_spvisc_ok) {
        updateSpeciesViscosities();
    }

    // see Eq. (9-5.15) of Reid, Prausnitz, and Poling
    int j, k;
    for (j = 0; j < m_nsp; j++) {
        for (k = j; k < m_nsp; k++) {
            vratiokj = m_visc[k]/m_visc[j];
            wratiojk = m_mw[j]/m_mw[k];

            // Note that m_wratjk(k,j) holds the square root of
            // m_wratjk(j,k)!
            factor1 = 1.0 + (m_sqvisc[k]/m_sqvisc[j]) * m_wratjk(k,j);
            m_phi(k,j) = factor1*factor1 /
                         (SqrtEight * m_wratkj1(j,k));
            m_phi(j,k) = m_phi(k,j)/(vratiokj * wratiojk);
        }
    }
    m_viscwt_ok = true;
}

}
