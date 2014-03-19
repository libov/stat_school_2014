// ROOT headers
#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TLegend.h>

// system headers
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <getopt.h>
using namespace std;

#include <TMassPeakFit.h>

Double_t gauss(Double_t *x, Double_t *par) {

    Double_t mean   = par[0];
    Double_t rms    = par[1];
    Double_t norm   = par[2];

    Double_t arg = x[0];

    return ( norm * exp( -0.5 * pow((arg -mean)/rms, 2) ) );
}

Double_t normalised_gauss(Double_t *x, Double_t *par) {

    Double_t mean   = par[0];
    Double_t rms    = par[1];
    Double_t area   = par[2];

    Double_t arg = x[0];

    Double_t prefactor = area / ( rms * sqrt(2*TMath::Pi()) );
    return ( prefactor * exp( -0.5 * pow((arg -mean)/rms, 2) ) );
}

Double_t pol3(Double_t *x, Double_t *par) {

    Double_t p0 = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];

    Double_t arg = x[0];

    return ( p0 + p1 * arg + p2 * pow(arg,2) + p3 * pow(arg,3) );
}

unsigned get_n_parameters(TString function) {

    map<TString,unsigned> n_par;

    n_par["gauss"] = 3;
    n_par["normalised_gauss"] = 3;
    n_par["pol3"]  = 4;

    return (n_par[function]);
}

typedef Double_t (*FIT_FUNCTION)(Double_t *x, Double_t *par);

FIT_FUNCTION get_function_pointer(TString function){ // another way: Double_t (*get_function_pointer(TString function))(Double_t *x, Double_t *par) 

    map<TString, FIT_FUNCTION> pointers;

    pointers["gauss"]   = gauss;
    pointers["normalised_gauss"]   = normalised_gauss;
    pointers["pol3"]    = pol3;

    return pointers[function];
}

TString get_parameter_name(TString function, unsigned par_nr) {

    map<TString, vector<TString> > parameter_names;

    parameter_names["gauss"].clear();
    parameter_names["gauss"].push_back("Mean             ");
    parameter_names["gauss"].push_back("RMS (sigma)      ");
    parameter_names["gauss"].push_back("Value at Maximum ");

    parameter_names["normalised_gauss"].clear();
    parameter_names["normalised_gauss"].push_back("Mean             ");
    parameter_names["normalised_gauss"].push_back("RMS (sigma)      ");
    parameter_names["normalised_gauss"].push_back("Area             ");

    parameter_names["pol3"].clear();
    parameter_names["pol3"].push_back("p0               ");
    parameter_names["pol3"].push_back("p1               ");
    parameter_names["pol3"].push_back("p2               ");
    parameter_names["pol3"].push_back("p3               ");

    return parameter_names[function][par_nr];
}
