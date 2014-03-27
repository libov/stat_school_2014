#ifndef TMASSPEAKFIT_H
#define TMASSPEAKFIT_H

#include <TROOT.h>
#include <TH1F.h>
#include <TMinuit.h>

#include <vector>
using namespace std;

class TMassPeakFit : public TObject {

    public:

        TMassPeakFit(TString config, void (*fcn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t));
        ~TMassPeakFit(){};

        void        CustomFit();
        void        Fit();
        Double_t    MinimizationFunction(Double_t * par);
        Double_t    FitFunction(Double_t * x, Double_t * par);
        void        MakePlots();
        Double_t    GetBinWidth() {return fBinWidth;}
        Double_t    GetHistogramLowerLimit() {return fHistogram -> GetBinLowEdge(1);}
        Double_t    GetHistogramUpperLimit() {return ( fHistogram -> GetBinLowEdge(fHistogram -> GetNbinsX()) + GetBinWidth() ) ;}
        bool        IsNormalisedGauss() {return fIsNormalisedGauss;}

    private:

        void        ReadSettings();

        // various minimization functions
        Double_t    chi2(Double_t * par);
        Double_t    log_likelihood(Double_t * par);

        // config file name
        TString     fConfigFilename;

        // data read off the config file
        TString     fFitFunction;
        unsigned    fNParameters;
        TString     fFitType;
        TString     fInputFile;
        TString     fInputHistogramName;
        Double_t    fStartValues[99];
        Double_t    fFitStep[99];
        Double_t    fLowerLimit[99];
        Double_t    fUpperLimit[99];
        Double_t    fFitRange[99];
        Double_t    fFixParameters[99];
        bool        fOnlyDrawInitial;

        unsigned        fNFitFunctions;
        vector<TString> fFitFunctions;
        bool            fPearson;
        bool            fDrawCovarianceEllipses;
        unsigned        fCovarianceEllipseParameter1;
        unsigned        fCovarianceEllipseParameter2;
        TString         fParName[99];
        double          fBinWidth;
        bool            fIsNormalisedGauss;

        // internal objects
        TH1F*       fHistogram;
        TMinuit*    fMinuit;

};

extern TMassPeakFit *gTMassPeakFit;

#endif
