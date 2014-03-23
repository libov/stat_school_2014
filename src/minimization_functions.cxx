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
        // check if we are in the fit range given in the config file
        if ( (x<fFitRange[0]) || (x>fFitRange[1]) ) continue;
        // calculate the term and add it to the overall chi2
        Double_t variance;
        if (fPearson) { // Pearson's chi2
            variance = FitFunction(&x, par);
        } else {  // Neyman's chi2
            variance = fHistogram -> GetBinContent(i);
            // sanity check
            if ( sqrt(variance) != fHistogram -> GetBinError(i) ) {
                cout << "WARNING: bin content error does not match sqrt(bin content). Should not happen unless the histogram was reweighted (was it?)." << endl;
                cout << "sqrt(bin content)= " << variance <<  ", bin error= " << fHistogram -> GetBinError(i) << endl;
            }
        }
        // this affects mostly Neyman's definition; entries with zero number of events should be skipped - the error is otherwise zero and can't be used in the denominator!
        if (variance == 0 ) continue;
        result += pow( fHistogram -> GetBinContent(i) - FitFunction(&x, par), 2 ) / variance;
    }

    // done
    return result;
}

// this is the Poisson's likelihood minimization function
Double_t TMassPeakFit::log_likelihood(Double_t * par) {

    // get number of bins in the input histogram - this gives number of terms in the chi2 calculation (hence the loop below)
    unsigned nbins = fHistogram -> GetNbinsX();

    // resulting likelihood will be stored here
    Double_t L = 1;

    // this is the value at which the fit function is evaluated
    Double_t x;

    // loop over factors entering L (=histogram bins)
    for (unsigned i=1; i<=nbins; i++) {
        // get the center of the bin
        x = fHistogram ->  GetBinCenter(i);
        // check if we are in the fit range given in the config file
        if ( (x<fFitRange[0]) || (x>fFitRange[1]) ) continue;
        // calculate the term and multiply the overall likelihood by it
        L *= exp(-FitFunction(&x, par)) * pow(FitFunction(&x, par), fHistogram -> GetBinContent(i)); // / TMath::Factorial(fHistogram -> GetBinContent(i));
    }

    // done
    return (-2*TMath::Log(L));
}
