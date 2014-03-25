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

    // if only_draw_initial is 1 in the config file, the minimization is not performed - instead function using start values are displayed
    if (fOnlyDrawInitial) return;

    // ***********  USER EDITABLE PART  *********** //
    
    // perform the minimization!
    fMinuit -> Command("MIGRAD");
    
    
    // ********  END OF USER EDITABLE PART  ******* //
}
