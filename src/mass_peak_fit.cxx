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

#include <TMassPeakFit.h>

void minimization_function(Int_t& npar, Double_t* grad, Double_t& f, Double_t par[], Int_t iflag) {
	f = gTMassPeakFit -> MinimizationFunction(par);
}

int main (int argc, char **argv) {

    // --------------------------------------------------------- //
    // --------------- read command line options --------------- //
    // --------------------------------------------------------- //

    static struct option long_options[] = {
        {"config",     required_argument, 0, 1},
        {"custom_fit",     no_argument,       0, 2},
    };
    TString  config     = "";
    bool     custom_fit = false;

    int option;
    int option_index;
    while ((option = getopt_long (argc, argv, "h", long_options, &option_index)) != -1) {
        switch (option) {
            case 1:
                config = optarg;
                cout << "\nINFO: config file is chosen to be ./config/mass_peak_fit/" << config << endl << endl;
                system("cat config/mass_peak_fit/" + config);
                cout << endl;
                break;
            case 2:
                custom_fit = true;
                cout << "INFO: working in custom mode. TMassPeakFit::CustomFit() will be executed" << endl;
                break;
            case 'h':
                cout<<"\nUsage: " << endl;
                cout<<"\tmass_peak_fit  --config <configuration file> [--custom_fit]\n"<<endl;
                cout << endl;
                cout << " AVAILABLE CONFIGURATION FILES:\n";
                system("ls config/mass_peak_fit | grep -v README");

                exit(-1);
            default:
                cout << "Unknown opiton or missing option argument. The program will terminate, sorry." << endl;
                exit(-1);
        }
    }

    // check if config file was set
    if (config=="") {
        cout << "ERROR: configuration file not set" << endl;
        cout << "Try\nmass_peak_fit -h\for more information" << endl;
        abort();
    }

    // create object
    TMassPeakFit instance(config, minimization_function);
    gTMassPeakFit = &instance;

    if (custom_fit) {
        instance.CustomFit();
    } else {
        instance.Fit();
    }

    instance.MakePlots();

    // done
    return 0;
}
