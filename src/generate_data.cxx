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
#include <TRandom3.h>
#include <TF1.h>

// system headers
#include <iostream>
#include<fstream>
#include <vector>
#include <getopt.h>
using namespace std;

struct Distribution {
    TString     type;
    unsigned    npoints;
    unsigned    seed;
    Double_t    par[99];
};

int main (int argc, char **argv) {

    // --------------------------------------------------------- //
    // --------------- read command line options --------------- //
    // --------------------------------------------------------- //

    static struct option long_options[] = {
        {"config",     required_argument, 0, 1},
    };
    TString  config     = "";

    int option;
    int option_index;
    while ((option = getopt_long (argc, argv, "h", long_options, &option_index)) != -1) {
        switch (option) {
            case 1:
                config = optarg;
                cout << "\nINFO: config file is chosen to be ./config/generate_data/" << config << endl << endl;
                system("cat config/generate_data/" + config);
                cout << endl << endl;
                break;
            case 'h':
                cout<<"\nUsage: " << endl;
                cout<<"\tgenerate_data --config <configuration file>\n\n\n";
                cout << " AVAILABLE CONFIGURATION FILES:\n";
                system("ls config/generate_data | grep -v README");
                exit(-1);
            default:
                cout << "Unknown opiton or missing option argument. The program will terminate, sorry." << endl;
                exit(-1);
        }
    }

    // check if config file was set
    if (config=="") {
        cout << "ERROR: configuration file not set" << endl;
        cout << "Try\ngenerate_data -h\for more information" << endl;
        abort();
    }

    // open config file
    TString config_path = "config/generate_data/"+config;
    ifstream f(config_path);
    if (!f.is_open()) {
        cout << "ERROR: Unable to open file " << config_path << endl;
        abort();
    }

    // --------------------------------------------------------- //
    // ------------ read settings from config file ------------- //
    // --------------------------------------------------------- //

    TString  output_file    = "";
    TString  histogram      = "";
    unsigned nbins          = 0;
    unsigned xmin           = 0;
    unsigned xmax           = 0;
    vector <Distribution> distributions;
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
        if (first_word == "output_file")    output_file = ((TObjString*)tokens->At(1)) -> GetString();
        if (first_word == "histogram")      histogram = ((TObjString*)tokens->At(1)) -> GetString();
        if (first_word == "nbins")          nbins = (((TObjString*)tokens->At(1)) -> GetString()).Atoi();
        if (first_word == "xmin")           xmin = (((TObjString*)tokens->At(1)) -> GetString()).Atoi();
        if (first_word == "xmax")           xmax = (((TObjString*)tokens->At(1)) -> GetString()).Atoi();
 
        if (first_word == "distribution") {
            Distribution distr;
            distr.type = ((TObjString*)tokens->At(1)) -> GetString();
            distr.npoints = (((TObjString*)tokens->At(2)) -> GetString()).Atoi();
            distr.seed = (((TObjString*)tokens->At(3)) -> GetString()).Atoi();
            unsigned entries = tokens -> GetEntries();
            unsigned npar = entries - 4; // total # of entries is # of function-specific parameters + 4 (keyword, function type, number of points, seed)
            for (unsigned i=0; i<npar; i++) {
                distr.par[i] = (((TObjString*)tokens->At(i+4)) -> GetString()).Atof();
            }
            distributions.push_back(distr);
        }
    }

    // --------------------------------------------------------- //
    // ------------- generate data, store to file-- ------------ //
    // --------------------------------------------------------- //

    // open output file
    TFile * output = new TFile ("data/"+output_file, "recreate");
    cout << "INFO: storing data to " << "data/"+output_file << endl;
    TH1F * h = new TH1F (histogram, "", nbins, xmin, xmax);
    
    for (unsigned i=0; i<distributions.size(); i++) {
        Distribution distr = distributions[i];
        TF1 *f;  // helping function for distributions not directly implemented in TRandom3
        Double_t normalisation;
        if (distr.type == "pol3") {
            f = new TF1("","pol3",xmin, xmax);
            f -> SetParameters(distr.par);
            normalisation = f -> Integral(xmin, xmax);
        }
        TRandom3 rnd(distr.seed);

        unsigned j=0;
        while (j<distr.npoints) {
            if (distr.type == "gauss") {
                double mean = distr.par[0];
                double width = distr.par[1];
                h -> Fill( rnd.Gaus(mean, width) );
                j++;
            } else if (distr.type == "pol3") {
                double x = rnd.Uniform(xmin, xmax);
                double y = rnd.Rndm();
                double max_y = f->Eval(x)/normalisation;
                if (y < max_y ) {
                    h -> Fill(x);
                    j++;
                }
            }
            
        }

        for (unsigned i=0; i<distr.npoints; i++) {
        }
    }

    h -> Write();
    output -> Close();

    // done
    cout << endl;
    return 0;
}
