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
    n_par["pol3"]  = 4;

    return (n_par[function]);
}

typedef Double_t (*FIT_FUNCTION)(Double_t *x, Double_t *par);

FIT_FUNCTION get_function_pointer(TString function){ // another way: Double_t (*get_function_pointer(TString function))(Double_t *x, Double_t *par) 

    map<TString, FIT_FUNCTION> pointers;

    pointers["gauss"]   = gauss;
    pointers["pol3"]    = pol3;

    return pointers[function];
}

TString get_parameter_name(TString function, unsigned par_nr) {

    TString result = "";

    if (function == "gauss") {

        if (par_nr == 0) result = "Mean            ";
        if (par_nr == 1) result = "RMS (sigma)     ";
        if (par_nr == 2) result = "Value at Maximum";

    } else if (function == "gauss+pol3") {

        if (par_nr == 0) result = "Mean            ";
        if (par_nr == 1) result = "RMS (sigma)     ";
        if (par_nr == 2) result = "Value at Maximum";
        if (par_nr == 3) result = "p0              ";
        if (par_nr == 4) result = "p1              ";
        if (par_nr == 5) result = "p2              ";
        if (par_nr == 6) result = "p3              ";

    } else if (function == "gauss+gauss+pol3") {

        if (par_nr == 0) result = "Mean (gauss1)            ";
        if (par_nr == 1) result = "RMS (sigma) (gauss1)     ";
        if (par_nr == 2) result = "Value at Maximum (gauss1)";
        if (par_nr == 3) result = "Mean (gauss2)            ";
        if (par_nr == 4) result = "RMS (sigma) (gauss2)     ";
        if (par_nr == 5) result = "Value at Maximum (gauss2)";
        if (par_nr == 6) result = "p0                       ";
        if (par_nr == 7) result = "p1                       ";
        if (par_nr == 8) result = "p2                       ";
        if (par_nr == 9) result = "p3                       ";

    } else {
        cout << "ERROR: fit function " << function << " is not supported " << endl;
    }

    return result;

}
