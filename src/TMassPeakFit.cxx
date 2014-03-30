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
#include <TPaveText.h>

// system headers
#include <iostream>
#include<fstream>
#include <vector>
#include <getopt.h>
using namespace std;

TMassPeakFit *gTMassPeakFit;

// ========== helping functions ========== //
TString toStr(Double_t arg, Int_t decimals);
TString toStr(Int_t arg);

TString toStr(Double_t arg, Int_t decimals) {
    char tmp[256];
    TString format = "%."+toStr(decimals)+"f";
    sprintf(tmp, format, arg);
    TString result(tmp);
    return result;
}

TString toStr(Int_t arg) {
    TString result;
    result += arg;
    return result;
}

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

    // initialise Minuit (set parameter names, start values, limits, fit steps; some of these are set in the config file, others directly here)
    for (unsigned i=0; i<fNParameters; i++) {
        fMinuit->DefineParameter(i, fParName[i], fStartValues[i], fFitStep[i], fLowerLimit[i], fUpperLimit[i]);
        if ( fFixParameters[i] == 1 ) fMinuit -> FixParameter(i);
    }
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
        if (first_word == "lower_limit") {
            for (Int_t i=0; i<tokens->GetEntries()-1; i++) {
                fLowerLimit[i] = (((TObjString*)tokens->At(i+1)) -> GetString()).Atof();
            }
        }
        if (first_word == "upper_limit") {
            for (Int_t i=0; i<tokens->GetEntries()-1; i++) {
                fUpperLimit[i] = (((TObjString*)tokens->At(i+1)) -> GetString()).Atof();
            }
        }
        if (first_word == "step") {
            for (Int_t i=0; i<tokens->GetEntries()-1; i++) {
                fFitStep[i] = (((TObjString*)tokens->At(i+1)) -> GetString()).Atof();
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
    fitresult -> SetNpx(1000);
    fitresult -> Draw("same");

    // draw individual components
    if (fNFitFunctions >1) {
        unsigned n_par_current=0;
        for (unsigned i=0; i<fNFitFunctions; i++) {
            TString function = fFitFunctions[i];
            TF1 * f = new TF1("", get_function_pointer(function), fFitRange[0], fFitRange[1], get_n_parameters(function));
            f -> SetParameters(&par[n_par_current]);
            n_par_current += get_n_parameters(function);
            f -> SetLineColor(i+3);
            f -> SetNpx(1000);
            f -> Draw("same");
        }
    }

    // and the legend
    TLegend * leg = new TLegend (0.15,0.6,0.35,0.8);
    leg -> AddEntry(fHistogram, "input data", "p");
    leg -> AddEntry(fitresult, "fit", "l");
    leg -> SetFillColor(0);
    leg -> Draw("same");
    
    // and fit results
    TPaveText *pt = new TPaveText(.7, .45, .98, .75, "NDC");
    pt -> SetFillColor(0);
    pt -> SetTextAlign(11);
    for (unsigned i=0; i<fNParameters; i++) {
        unsigned ndigits = 2;
        if ( fParName[i].Contains("Sigma") || fParName[i].Contains("Mu")) ndigits = 5;
        TString text = fParName[i] + " = " + toStr(par[i], ndigits) + " #pm " + toStr(par_err[i], ndigits);
        if ( fFixParameters[i] != 1 ) pt -> AddText(text);
    }

    // print chi2
    bool tmp = fPearson;
    fPearson = false;
    double chi2_neyman = chi2(par);
    fPearson = tmp;
    unsigned ndof = fNDataPoints - gMinuit ->  GetNumFreePars();
    pt -> AddText("#chi^{2}/ndof = " + toStr(chi2_neyman, 1) + "/" + toStr(ndof) + " = " + toStr(chi2_neyman/ndof, 2) );
    pt -> Draw("same");

    // print to file
    c -> Print("results/fit.eps");
    c -> Print("results/fit.root");

    // covariance ellipse
    if (fDrawCovarianceEllipses) {

        TCanvas * c2 = new TCanvas();
 
        gMinuit -> SetErrorDef(3*3);
        TGraph *cont1 = (TGraph*)gMinuit->Contour(100, fCovarianceEllipseParameter1-1, fCovarianceEllipseParameter2-1);
        if (cont1) cont1->Draw("al");

        gMinuit -> SetErrorDef(2.5*2.5);
        TGraph *cont2 = (TGraph*)gMinuit->Contour(100, fCovarianceEllipseParameter1-1, fCovarianceEllipseParameter2-1);
        if (cont2) cont2->Draw("samel");
        
        gMinuit -> SetErrorDef(2*2);
        TGraph *cont3 = (TGraph*)gMinuit->Contour(100, fCovarianceEllipseParameter1-1, fCovarianceEllipseParameter2-1);
        if (cont3) cont3->Draw("samel");
        
        gMinuit -> SetErrorDef(1.5*1.5);
        TGraph *cont4 = (TGraph*)gMinuit->Contour(100, fCovarianceEllipseParameter1-1, fCovarianceEllipseParameter2-1);
        if (cont4) cont4->Draw("samel");

        gMinuit -> SetErrorDef(1*1);
        TGraph *cont5 = (TGraph*)gMinuit->Contour(100, fCovarianceEllipseParameter1-1, fCovarianceEllipseParameter2-1);
        if (cont5) cont5->Draw("samel");
        
        gMinuit -> SetErrorDef(0.5*0.5);
        TGraph *cont6 = (TGraph*)gMinuit->Contour(100, fCovarianceEllipseParameter1-1, fCovarianceEllipseParameter2-1);
        if (cont6) cont6->Draw("samel");

        if (cont1) {
  
            cont1 -> SetTitle ("Covariance ellipse");
            cont1 -> GetXaxis() -> SetTitle(fParName[fCovarianceEllipseParameter1-1]);
            cont1 -> GetYaxis() -> SetTitle(fParName[fCovarianceEllipseParameter2-1]);
            cont1 -> GetXaxis() -> SetTitleOffset(1.5);
            cont1 -> GetYaxis() -> SetTitleOffset(1.5);

            c2 -> Print("results/cov_ellipse.eps");
            c2 -> Print("results/cov_ellipse.root");
        }
    }

    // pulls
    TH1F * pull1 = new TH1F ("", "", 20, -3, 3);
    double x_down = GetHistogramLowerLimit();
    double x_up = GetHistogramUpperLimit();
    unsigned nbins = fHistogram -> GetNbinsX();
    TH1F * pull2 = new TH1F ("", "", nbins, x_down, x_up);
    for (unsigned i=1; i<=nbins; i++) {
        double data = fHistogram -> GetBinContent(i);
        double fit = fitresult -> Eval(fHistogram -> GetBinCenter(i));
        double deviation = data - fit;
        double error = fHistogram -> GetBinError(i);
        if (data!=0) {
            pull1 -> Fill(deviation/error);
            pull2 -> SetBinContent(i, deviation);
            pull2 -> SetBinError(i, error);
        }
    }

    TCanvas c3;
    c3.cd();
    pull1 -> Draw();
    c3.Print("results/pull1.eps");
    c3.Print("results/pull1.root");

    TCanvas c4;
    c4.cd();
    pull2 -> Draw();
    TLine * l = new TLine(x_down, 0, x_up, 0);
    l -> SetLineStyle(2);
    l -> Draw("same");
    c4.Print("results/pull2.eps");
    c4.Print("results/pull2.root");
}
