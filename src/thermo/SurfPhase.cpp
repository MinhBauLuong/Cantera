/**
 *  @file SurfPhase.cpp
 *  Definitions for a simple thermodynamic model of a surface phase
 *  derived from ThermoPhase,  assuming an ideal solution model
 *  (see \ref thermoprops and class
 *  \link Cantera::SurfPhase SurfPhase\endlink).
 */
// Copyright 2002  California Institute of Technology

#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/EdgePhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

using namespace ctml;
using namespace std;

namespace Cantera
{
SurfPhase::SurfPhase(doublereal n0):
    ThermoPhase(),
    m_n0(n0),
    m_logn0(0.0),
    m_press(OneAtm),
    m_tlast(0.0)
{
    if (n0 > 0.0) {
        m_logn0 = log(n0);
    }
    setNDim(2);
}

SurfPhase::SurfPhase(const std::string& infile, std::string id_) :
    ThermoPhase(),
    m_n0(0.0),
    m_logn0(0.0),
    m_press(OneAtm),
    m_tlast(0.0)
{
    XML_Node* root = get_XML_File(infile);
    if (id_ == "-") {
        id_ = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id_, root);
    if (!xphase) {
        throw CanteraError("SurfPhase::SurfPhase",
                           "Couldn't find phase name in file:" + id_);
    }
    // Check the model name to ensure we have compatibility
    const XML_Node& th = xphase->child("thermo");
    string model = th["model"];
    if (model != "Surface" && model != "Edge") {
        throw CanteraError("SurfPhase::SurfPhase",
                           "thermo model attribute must be Surface or Edge");
    }
    importPhase(*xphase, this);
}

SurfPhase::SurfPhase(XML_Node& xmlphase) :
    ThermoPhase(),
    m_n0(0.0),
    m_logn0(0.0),
    m_press(OneAtm),
    m_tlast(0.0)
{
    const XML_Node& th = xmlphase.child("thermo");
    string model = th["model"];
    if (model != "Surface" && model != "Edge") {
        throw CanteraError("SurfPhase::SurfPhase",
                           "thermo model attribute must be Surface or Edge");
    }
    importPhase(xmlphase, this);
}

SurfPhase::SurfPhase(const SurfPhase& right) :
    m_n0(right.m_n0),
    m_logn0(right.m_logn0),
    m_press(right.m_press),
    m_tlast(right.m_tlast)
{
    *this = operator=(right);
}

SurfPhase& SurfPhase::
operator=(const SurfPhase& right)
{
    if (&right != this) {
        ThermoPhase::operator=(right);
        m_n0         = right.m_n0;
        m_logn0      = right.m_logn0;
        m_press      = right.m_press;
        m_tlast      = right.m_tlast;
        m_h0         = right.m_h0;
        m_s0         = right.m_s0;
        m_cp0        = right.m_cp0;
        m_mu0        = right.m_mu0;
        m_work       = right.m_work;
        m_logsize    = right.m_logsize;
    }
    return *this;
}

ThermoPhase* SurfPhase::duplMyselfAsThermoPhase() const
{
    return new SurfPhase(*this);
}

doublereal SurfPhase::enthalpy_mole() const
{
    if (m_n0 <= 0.0) {
        return 0.0;
    }
    _updateThermo();
    return mean_X(DATA_PTR(m_h0));
}

doublereal SurfPhase::intEnergy_mole() const
{
    return enthalpy_mole();
}

void SurfPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    doublereal rt = GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= rt;
    }
}

void SurfPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R(sbar);
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }
}

void SurfPhase::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

// HKM 9/1/11  The partial molar volumes returned here are really partial molar areas.
//             Partial molar volumes for this phase should actually be equal to zero.
void SurfPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

void SurfPhase::getStandardChemPotentials(doublereal* mu0) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), mu0);
}

void SurfPhase::getChemPotentials(doublereal* mu) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), mu);
    getActivityConcentrations(DATA_PTR(m_work));
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += GasConstant * temperature() *
                 (log(m_work[k]) - logStandardConc(k));
    }
}

void SurfPhase::getActivityConcentrations(doublereal* c) const
{
    getConcentrations(c);
}

doublereal SurfPhase::standardConcentration(size_t k) const
{
    return m_n0/size(k);
}

doublereal SurfPhase::logStandardConc(size_t k) const
{
    return m_logn0 - m_logsize[k];
}

void SurfPhase::setParameters(int n, doublereal* const c)
{
    warn_deprecated("SurfPhase::setParameters");
    if (n != 1) {
        throw CanteraError("SurfPhase::setParameters",
                           "Bad value for number of parameter");
    }
    setSiteDensity(c[0]);
}

void SurfPhase::getGibbs_RT(doublereal* grt) const
{
    _updateThermo();
    double rrt = 1.0/(GasConstant*temperature());
    scale(m_mu0.begin(), m_mu0.end(), grt, rrt);
}

void SurfPhase::
getEnthalpy_RT(doublereal* hrt) const
{
    _updateThermo();
    double rrt = 1.0/(GasConstant*temperature());
    scale(m_h0.begin(), m_h0.end(), hrt, rrt);
}

void SurfPhase::getEntropy_R(doublereal* sr) const
{
    _updateThermo();
    double rr = 1.0/GasConstant;
    scale(m_s0.begin(), m_s0.end(), sr, rr);
}

void SurfPhase::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    double rr = 1.0/GasConstant;
    scale(m_cp0.begin(), m_cp0.end(), cpr, rr);
}

void SurfPhase::getStandardVolumes(doublereal* vol) const
{
    _updateThermo();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = 1.0/standardConcentration(k);
    }
}

void SurfPhase::getGibbs_RT_ref(doublereal* grt) const
{
    getGibbs_RT(grt);
}

void SurfPhase::getEnthalpy_RT_ref(doublereal* hrt) const
{
    getEnthalpy_RT(hrt);
}

void SurfPhase::getEntropy_R_ref(doublereal* sr) const
{
    getEntropy_R(sr);
}

void SurfPhase::getCp_R_ref(doublereal* cprt) const
{
    getCp_R(cprt);
}

void SurfPhase::initThermo()
{
    if (m_kk == 0) {
        throw CanteraError("SurfPhase::initThermo",
                           "Number of species is equal to zero");
    }
    m_h0.resize(m_kk);
    m_s0.resize(m_kk);
    m_cp0.resize(m_kk);
    m_mu0.resize(m_kk);
    m_work.resize(m_kk);
    vector_fp cov(m_kk, 0.0);
    cov[0] = 1.0;
    setCoverages(DATA_PTR(cov));
    m_logsize.resize(m_kk);
    for (size_t k = 0; k < m_kk; k++) {
        m_logsize[k] = log(size(k));
    }
}

void SurfPhase::setSiteDensity(doublereal n0)
{
    if (n0 <= 0.0) {
        throw CanteraError("SurfPhase::setSiteDensity",
                           "Bad value for parameter");
    }
    m_n0 = n0;
    m_logn0 = log(m_n0);
}

void SurfPhase::setCoverages(const doublereal* theta)
{
    double sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += theta[k];
    }
    if (sum <= 0.0) {
        for (size_t k = 0; k < m_kk; k++) {
            cout << "theta(" << k << ") = " << theta[k] << endl;
        }
        throw CanteraError("SurfPhase::setCoverages",
                           "Sum of Coverage fractions is zero or negative");
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_work[k] = m_n0*theta[k]/(sum*size(k));
    }
    /*
     * Call the Phase:: class function
     * setConcentrations.
     */
    setConcentrations(DATA_PTR(m_work));
}

void SurfPhase::setCoveragesNoNorm(const doublereal* theta)
{
    for (size_t k = 0; k < m_kk; k++) {
        m_work[k] = m_n0*theta[k]/(size(k));
    }
    /*
     * Call the Phase:: class function
     * setConcentrations.
     */
    setConcentrations(DATA_PTR(m_work));
}

void SurfPhase::getCoverages(doublereal* theta) const
{
    getConcentrations(theta);
    for (size_t k = 0; k < m_kk; k++) {
        theta[k] *= size(k)/m_n0;
    }
}

void SurfPhase::setCoveragesByName(const std::string& cov)
{
    size_t kk = nSpecies();
    compositionMap cc = parseCompString(cov, speciesNames());
    doublereal c;
    vector_fp cv(kk, 0.0);
    bool ifound = false;
    for (size_t k = 0; k < kk; k++) {
        c = cc[speciesName(k)];
        if (c > 0.0) {
            ifound = true;
            cv[k] = c;
        }
    }
    if (!ifound) {
        throw CanteraError("SurfPhase::setCoveragesByName",
                           "Input coverages are all zero or negative");
    }
    setCoverages(DATA_PTR(cv));
}

void SurfPhase::_updateThermo(bool force) const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow || force) {
        m_spthermo->update(tnow, DATA_PTR(m_cp0), DATA_PTR(m_h0),
                           DATA_PTR(m_s0));
        m_tlast = tnow;
        doublereal rt = GasConstant * tnow;
        for (size_t k = 0; k < m_kk; k++) {
            m_h0[k] *= rt;
            m_s0[k] *= GasConstant;
            m_cp0[k] *= GasConstant;
            m_mu0[k] = m_h0[k] - tnow*m_s0[k];
        }
        m_tlast = tnow;
    }
}

void SurfPhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","Surface");
    doublereal n = getFloat(eosdata, "site_density", "toSI");
    if (n <= 0.0)
        throw CanteraError("SurfPhase::setParametersFromXML",
                           "missing or negative site density");
    m_n0 = n;
    m_logn0 = log(m_n0);
}

void SurfPhase::setStateFromXML(const XML_Node& state)
{

    double t;
    if (getOptionalFloat(state, "temperature", t, "temperature")) {
        setTemperature(t);
    }

    if (state.hasChild("coverages")) {
        string comp = getChildValue(state,"coverages");
        setCoveragesByName(comp);
    }
}

EdgePhase::EdgePhase(doublereal n0) : SurfPhase(n0)
{
    setNDim(1);
}

EdgePhase::EdgePhase(const EdgePhase& right) :
    SurfPhase(right.m_n0)
{
    setNDim(1);
    *this = operator=(right);
}

EdgePhase& EdgePhase::operator=(const EdgePhase& right)
{
    if (&right != this) {
        SurfPhase::operator=(right);
        setNDim(1);
    }
    return *this;
}

ThermoPhase* EdgePhase::duplMyselfAsThermoPhase() const
{
    return new EdgePhase(*this);
}

void EdgePhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","Edge");
    doublereal n = getFloat(eosdata, "site_density", "toSI");
    if (n <= 0.0)
        throw CanteraError("EdgePhase::setParametersFromXML",
                           "missing or negative site density");
    m_n0 = n;
    m_logn0 = log(m_n0);
}

}
