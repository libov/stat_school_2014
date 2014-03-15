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

    private:

        void        ReadSettings();

        // various minimization functions
        Double_t    chi2(Double_t * par);
        Double_t    poisson(Double_t * par);

        // config file name
        TString     fConfigFilename;

        // data read off the config file
        TString     fFitFunction;
        unsigned    fNParameters;
        TString     fFitType;
        TString     fInputFile;
        TString     fInputHistogramName;
        Double_t    fStartValues[99];
        Double_t    fFitRange[99];
        Double_t    fFixParameters[99];
        bool        fOnlyDrawInitial;

        // internal objects
        TH1F*       fHistogram;
        TMinuit*    fMinuit;

};
#endif
