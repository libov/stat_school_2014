#ifndef FITFUNCTIONS_H
#define FITFUNCTIONS_H

#include<TROOT.h>

Double_t gauss(Double_t *x, Double_t *par);
Double_t breit_wigner(Double_t *x, Double_t *par);

Double_t pol3(Double_t *x, Double_t *par);

TString get_parameter_name(TString function, unsigned par_nr);

#endif
