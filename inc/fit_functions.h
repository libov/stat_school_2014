#ifndef FITFUNCTIONS_H
#define FITFUNCTIONS_H

#include<TROOT.h>

Double_t gauss(Double_t *x, Double_t *par);

Double_t pol3(Double_t *x, Double_t *par);

TString get_parameter_name(TString function, unsigned par_nr);

unsigned get_n_parameters(TString function);

typedef Double_t (*FIT_FUNCTION)(Double_t *x, Double_t *par);

FIT_FUNCTION get_function_pointer(TString function); // another way: Double_t (*get_function_pointer(TString function))(Double_t *x, Double_t *par);

#endif
