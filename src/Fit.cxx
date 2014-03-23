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
#include <TGraph.h>

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

    // initialise Minuit (set parameter names, start values, limits, fit steps; some of these are set in the config file, others directly here)
    Double_t fitStep[99], limitMin[99], limitMax[99];
    for (int i=0; i<99; i++) {
        fitStep[i]  = 0.01;
        limitMin[i] = 0;
        limitMax[i] = 0;
    }

    for (unsigned i=0;i<fNParameters;i++) {
        fMinuit->DefineParameter(i, fParName[i], fStartValues[i], fitStep[i], limitMin[i], limitMax[i]);
    }

    // if only_draw_initial is 1 in the config file, the minimization is not performed - instead function using start values are displayed
    if (fOnlyDrawInitial) return;

    // perform the minimization!
    fMinuit -> Migrad();
}
