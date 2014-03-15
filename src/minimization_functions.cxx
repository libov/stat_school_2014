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

// Minimization functions for Minuit are implemented in this file
// Function names correspond to the option in the config file
// fHistogram is pointer to the histogram containing input data

// This is the chi2 minimization function
Double_t TMassPeakFit::chi2(Double_t * par) {

    // get number of bins in the input histogram - this gives number of terms in the chi2 calculation (hence the loop below)
    unsigned nbins = fHistogram -> GetNbinsX();

    // resulting chi2 will be stored here
    Double_t result = 0;

    // this is the value at which the fit function is evaluated
    Double_t x;

    // loop over chi2 terms (=histogram bins)
    for (unsigned i=1; i<=nbins; i++) {
        // get the center of the bin
        x = fHistogram ->  GetBinCenter(i);
        // skip for zero entries bins - the error is otherwise zero!
        if ( fHistogram -> GetBinContent(i) == 0 ) continue;
        // check if we are in the fit range given in the config file
        if ( (x<fFitRange[0]) || (x>fFitRange[1]) ) continue;
        // calculate the term and add it to the overall chi2
        result += pow( (fHistogram -> GetBinContent(i) - FitFunction(&x, par) )/fHistogram->GetBinError(i), 2);
    }

    // done
    return result;
}
