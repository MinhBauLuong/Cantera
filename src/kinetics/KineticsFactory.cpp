/**
 *  @file KineticsFactory.cpp
 */
// Copyright 2001  California Institute of Technology

#include "cantera/kinetics/KineticsFactory.h"

#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/GRI_30_Kinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/EdgeKinetics.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/kinetics/AqueousKinetics.h"

#include "cantera/kinetics/METHANE19Kinetics.h"
#include "cantera/kinetics/ETHANOL28Kinetics.h"
#include "cantera/kinetics/METHANE13Kinetics.h"
#include "cantera/kinetics/DME30Kinetics.h"
#include "cantera/kinetics/HEPTANE52Kinetics.h"
#include "cantera/kinetics/ISOOCTANE99Kinetics.h"
#include "cantera/kinetics/PRF116Kinetics.h"
#include  <iostream>

using namespace std;

namespace Cantera
{

KineticsFactory* KineticsFactory::s_factory = 0;
mutex_t KineticsFactory::kinetics_mutex;

static int ntypes = 13;
static string _types[] = {"none", "GasKinetics", "GRI30", "Interface", "Edge", "AqueousKinetics", "METHANE19", \
                          "ETHANOL28", "METHANE13", "DME30", "HEPTANE52", "ISOOCTANE99", "PRF116"};
static int _itypes[]   = {0, cGasKinetics, cGRI30, cInterfaceKinetics, cEdgeKinetics, cAqueousKinetics, cMETHANE19Kinetics,\
                          cETHANOL28Kinetics, cMETHANE13Kinetics, cDME30Kinetics, cHEPTANE52Kinetics, cISOOCTANE99Kinetics, cPRF116Kinetics};

Kinetics* KineticsFactory::
newKinetics(XML_Node& phaseData, vector<ThermoPhase*> th)
{
    /*
     * Look for a child of the xml element phase called
     * "kinetics". It has an attribute name "model".
     * Store the value of that attribute in the variable kintype
     */
    string kintype = phaseData.child("kinetics")["model"];
    /*
     * look up the string kintype in the list of known
     * kinetics managers (list is kept at the top of this file).
     * Translate it to an integer value, ikin.
     */
    int ikin=-1;
    int n;
    for (n = 0; n < ntypes; n++) {
        if (kintype == _types[n]) {
            ikin = _itypes[n];
        }
    }
    /*
     * Assign the kinetics manager based on the value of ikin.
     * Kinetics managers are classes derived from the base
     * Kinetics class. Unknown kinetics managers will throw a
     * CanteraError here.
     */
    Kinetics* k=0;
    switch (ikin) {

    case 0:
        k = new Kinetics;
        break;

    case cGasKinetics:
        k = new GasKinetics;
        break;

    case cGRI30:
        k = new GRI_30_Kinetics;
        break;

    case cInterfaceKinetics:
        k = new InterfaceKinetics;
        break;

    case cEdgeKinetics:
        k = new EdgeKinetics;
        break;

    case cAqueousKinetics:
        k = new AqueousKinetics;
        break;

    case cMETHANE19Kinetics:
        k = new METHANE19Kinetics;
        //printf("\nYou are using the Lu and Law 19 species analytical mechanism for CH4\n");
        break;

    case cETHANOL28Kinetics:
        k = new ETHANOL28Kinetics;
        //printf("You are blabla");
        break;

    case cMETHANE13Kinetics:
        k = new METHANE13Kinetics;
        //printf("\nYou are using the Lu and Law 19 species analytical mechanism for CH4\n");
        break;

    case cDME30Kinetics:
        k = new DME30Kinetics;
        break;

    case cHEPTANE52Kinetics:
        k = new HEPTANE52Kinetics;
        break;
    
    case cISOOCTANE99Kinetics:
        k = new ISOOCTANE99Kinetics;
        break;
    
    case cPRF116Kinetics:
        k = new PRF116Kinetics;
        break;
        
    default:
        throw UnknownKineticsModel("KineticsFactory::newKinetics",
                                   kintype);
    }

    // Now that we have the kinetics manager, we can
    // import the reaction mechanism into it.
    importKinetics(phaseData, th, k);

    // Return the pointer to the kinetics manager
    return k;
}

Kinetics* KineticsFactory::newKinetics(const string& model)
{

    int ikin = -1;
    int n;
    for (n = 0; n < ntypes; n++) {
        if (model == _types[n]) {
            ikin = _itypes[n];
        }
    }
    Kinetics* k=0;
    switch (ikin) {

    case cGasKinetics:
        k = new GasKinetics;
        break;

    case cGRI30:
        k = new GRI_30_Kinetics;
        break;

    case cInterfaceKinetics:
        k = new InterfaceKinetics;
        break;

    case cMETHANE19Kinetics:
        k = new METHANE19Kinetics;
        //printf("\nYou are using the Lu and Law 19 species analytical mechanism for CH4\n");
        break;

    case cETHANOL28Kinetics:
       k = new ETHANOL28Kinetics;
        //printf("You are blabla");
        break;

    case cMETHANE13Kinetics:
        k = new METHANE13Kinetics;
        //printf("\nYou are using the Lu and Law 19 species analytical mechanism for CH4\n");
        break;

    case cDME30Kinetics:
        k = new DME30Kinetics;
        //printf("\nYou are using the Lu and Law 19 species analytical mechanism for CH4\n");
        break;

    case cHEPTANE52Kinetics:
        k = new HEPTANE52Kinetics;
        //printf("\nYou are using the Lu and Law 19 species analytical mechanism for CH4\n");
        break;
      
    case cISOOCTANE99Kinetics:
        k = new ISOOCTANE99Kinetics;
        break;
    
    case cPRF116Kinetics:
        k = new PRF116Kinetics;
        break;
        
    default:
        throw UnknownKineticsModel("KineticsFactory::newKinetics",
                                   model);
    }
    return k;
}

}
