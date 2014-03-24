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

TMassPeakFit *gTMassPeakFit;

// --------------------------------------------------------- //
// ------------------ Class constructor -------------------- //
// --------------------------------------------------------- //
TMassPeakFit::TMassPeakFit(TString config, void (*fcn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t))
{
    // set the global pointer right away
    gTMassPeakFit = this;

    fDrawCovarianceEllipses = false;

    fConfigFilename = config;
    ReadSettings();

    // open input data
    TFile * input = new TFile ("data/"+fInputFile, "read");
    fHistogram = (TH1F*) input -> Get(fInputHistogramName);
    fBinWidth = fHistogram -> GetBinWidth(1);
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
	if (first_word == "normalised_gauss")   fIsNormalisedGauss = (bool) (((TObjString*)tokens->At(1)) -> GetString()).Atoi();
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
            fDrawCovarianceEllipses = true;
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

    unsigned k=0;
    for (unsigned i=0; i<fNFitFunctions; i++) {
        TString function = fFitFunctions[i];
        unsigned npar = get_n_parameters(function);
        for (unsigned j=0; j<npar; j++) {
            fParName[k] = get_parameter_name(function, j);
            k++;
        }
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

// --------------------------------------------------------- //
// ------------ produce plots using fit results ------------ //
// --------------------------------------------------------- //
void    TMassPeakFit::MakePlots() {
 
    TCanvas * c = new TCanvas();
    fHistogram -> Draw("e");
    fHistogram -> SetAxisRange(0, 1.2*fHistogram->GetMaximum(), "Y");

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
    c -> Print("results/fit.eps");
    c -> Print("results/fit.root");

    // draw covariance ellipse
    if (fDrawCovarianceEllipses) {
        TCanvas * c2 = new TCanvas();
        TGraph *cont = (TGraph*)gMinuit->Contour(100, fCovarianceEllipseParameter1, fCovarianceEllipseParameter2);
        if (cont) {
            cont -> SetTitle ("Covariance ellipse");
            cont -> GetXaxis() -> SetTitle(fParName[fCovarianceEllipseParameter1]);
            cont -> GetYaxis() -> SetTitle(fParName[fCovarianceEllipseParameter2]);
            cont->Draw("al");
            c2 -> Print("results/cov_ellipse.eps");
            c2 -> Print("results/cov_ellipse.root");
        }
    }
}
