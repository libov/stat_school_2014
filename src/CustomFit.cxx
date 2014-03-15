#include <TMassPeakFit.h>

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

// fill this function with your own stuff
// remember to use --custom_fit command line option with mass_peak_fit program to invoke this function!
void TMassPeakFit::CustomFit() {
    // TMinuit * Minuit = new TMinuit(3);
    // Minuit -> Migrad();
};
