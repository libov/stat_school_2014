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
// ------------------ Class constructor -------------------- //
// --------------------------------------------------------- //
TMassPeakFit::TMassPeakFit(TString config, void (*fcn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t))
{
    fConfigFilename = config;
    ReadSettings();

    // open input data
    TFile * input = new TFile ("data/"+fInputFile, "read");
    fHistogram = (TH1F*) input -> Get(fInputHistogramName);
    fHistogram -> SetMarkerStyle(20);

    fMinuit = new TMinuit(3);
    fMinuit -> SetFCN(fcn);
}

// --------------------------------------------------------- //
// ------------ read settings from config file ------------- //
// --------------------------------------------------------- //
void    TMassPeakFit::ReadSettings() {

    // open config file
    TString config_path = "config/mass_peak_fit/"+fConfigFilename;
    ifstream f(config_path);
    if (!f.is_open()) {
        cout << "ERROR: Unable to open file " << config_path << endl;
        abort();
    }
    cout << "\nINFO: opened " << config_path << endl;

    string line;
    while ( f.good() ) {

        // read each line
        getline (f,line);

        // skip if an empty line
        if (line=="") continue;
        // tokenize
        TString line_str = line;
        TObjArray * tokens = line_str.Tokenize(" ");

        // check if this line is a comment
        TString first_word = ((TObjString*)tokens->At(0)) -> GetString();
        char first_char = first_word[0];
        if (first_char=='#') continue;

        // depending on the keyword, store values to variables
        if (first_word == "fit_function")       fFitFunction = ((TObjString*)tokens->At(1)) -> GetString();
        if (first_word == "fit_type")           fFitType = ((TObjString*)tokens->At(1)) -> GetString();
        if (first_word == "input_file")         fInputFile = ((TObjString*)tokens->At(1)) -> GetString();
        if (first_word == "input_histogram")    fInputHistogramName = ((TObjString*)tokens->At(1)) -> GetString();
        if (first_word == "fit_range") {
            fFitRange[0] = (((TObjString*)tokens->At(1)) -> GetString()).Atof();
            fFitRange[1] = (((TObjString*)tokens->At(2)) -> GetString()).Atof();
        }
        if (first_word == "start_values") {
            fNParameters = tokens->GetEntries()-1;
            for (unsigned i=0; i<fNParameters; i++) {
                fStartValues[i] = (((TObjString*)tokens->At(i+1)) -> GetString()).Atof();
            }
        }
        if (first_word == "fix_parameters") {
            for (Int_t i=0; i<tokens->GetEntries()-1; i++) {
                fFixParameters[i] = (((TObjString*)tokens->At(i+1)) -> GetString()).Atof();
            }
        }
        if (first_word == "covariance_ellipse_parameters") {
            fCovarianceEllipseParameter1 =(((TObjString*)tokens->At(1)) -> GetString()).Atoi();
            fCovarianceEllipseParameter2 =(((TObjString*)tokens->At(2)) -> GetString()).Atoi();
        }
        if (first_word == "only_draw_initial") {
            fOnlyDrawInitial = (bool) (((TObjString*)tokens->At(1)) -> GetString()).Atoi();
        }
    }

    // store names of individual fit functions (which are summed up to form the overall fit function)
    TObjArray *tokens = fFitFunction.Tokenize("+");
    fNFitFunctions = tokens -> GetEntries();
    fFitFunctions.clear();
    unsigned n_parameters = 0;
    for (unsigned i=0; i<fNFitFunctions; i++) {
        TString function = ( (TObjString*)tokens->At(i) ) -> GetString();
         fFitFunctions.push_back(function);
         n_parameters += get_n_parameters(function);
    }

    // sanity check - number of parameters
    if (fNParameters != n_parameters ) {
        cout << "ERROR: function type and number of parameters given in start_values don't match ... " << endl;
        abort();
    }

    // this is to distinguish between Neyman's and Pearson's chi2 definitions
    fPearson = false;
    if (fFitType == "chi2_pearson") {
        fPearson = true;
    }
}

// --------------------------------------------------------- //
// ----------- Minimization Function for minuit ------------ //
// --------------------------------------------------------- //
Double_t    TMassPeakFit::MinimizationFunction(Double_t * par) {
    Double_t result;
    if ( (fFitType == "chi2_neyman") || (fFitType == "chi2_pearson") ) {
        result = chi2(par);
    } else if ( fFitType == "log_likelihood") {
        return log_likelihood(par);
    } else {
        cout << "ERROR: minimisation function " << fFitType << " is not supported" << endl;
        abort();
    }
    return result;
}

// --------------------------------------------------------- //
// --------------------- Fit Function  --------------------- //
// --------------------------------------------------------- //
Double_t    TMassPeakFit::FitFunction(Double_t * x, Double_t * par) {

    Double_t result = 0;

    unsigned npar = 0;

    for (unsigned i=0; i<fNFitFunctions; i++) {

        TString function = fFitFunctions[i];

        result += get_function_pointer(function)(x, &par[npar]);

        npar += get_n_parameters(function);
    }

    return result;
}
