#include <TMassPeakFit.h>
#include <fit_functions.h>

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
#include<fstream>
#include <vector>
#include <getopt.h>
using namespace std;

// --------------------------------------------------------- //
// ------------------- perform the fit --------------------- //
// --------------------------------------------------------- //
void TMassPeakFit::Fit() {

    TCanvas * c = new TCanvas();
    fHistogram -> Draw("e");
    fHistogram -> SetAxisRange(0, 1.2*fHistogram->GetMaximum(), "Y");

    // initialise Minuit (set parameter names, start values, limits, fit steps; some of these are set in the config file, others directly here.
    TString  ParName[99];
    for (unsigned i=0; i<fNParameters; i++) ParName[i] = get_parameter_name(fFitFunction, i);

    Double_t fitStep[99], limitMin[99], limitMax[99];

    for (int i=0; i<99; i++) {
        fitStep[i]  = 0.01;
        limitMin[i] = 0;
        limitMax[i] = 0;
    }

    for (unsigned i=0;i<fNParameters;i++) {
        fMinuit->DefineParameter(i, ParName[i], fStartValues[i], fitStep[i], limitMin[i], limitMax[i]);
    }

    // perform the minimization!
    if (!fOnlyDrawInitial) fMinuit -> Migrad();

    // get parameter values and their uncertainties from Minuit
    Double_t par[99];
    Double_t par_err[99];
    for (unsigned i=0; i<fNParameters; i++) {
        fMinuit -> GetParameter(i, par[i], par_err[i]);
    }

    // draw the fitted function
    TF1 * fitresult =  new TF1("f", this, &TMassPeakFit::FitFunction, fFitRange[0], fFitRange[1], fNParameters, "TMassPeakFit", "FitFunction");
    fitresult -> SetParameters(par);
    fitresult -> Draw("same");

    // draw individual components
    unsigned n_par_current=0;
    for (unsigned i=0; i<fNFitFunctions; i++) {
        TString function = fFitFunctions[i];
        TF1 * f = new TF1("", get_function_pointer(function), fFitRange[0], fFitRange[1], get_n_parameters(function));
        f -> SetParameters(&par[n_par_current]);
        n_par_current += get_n_parameters(function);
        f -> SetLineColor(i+3);
        f -> Draw("same");
    }

    // and the legend
    TLegend * leg = new TLegend (0.15,0.6,0.35,0.8);
    leg -> AddEntry(fHistogram, "input data", "p");
    leg -> AddEntry(fitresult, "fit", "l");
    leg -> SetFillColor(0);
    leg -> Draw("same");

    // print to file
    c -> Print("results/out.eps");
}
