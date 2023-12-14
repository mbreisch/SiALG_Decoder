#ifndef FITMODELS_H
#define FITMODELS_H

#include <SerialisableObject.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>

using namespace std;

class FitModels{

    friend class boost::serialization::access;

    public:

        FitModels();
        ~FitModels();

        static Double_t Gaussian(Double_t *x, Double_t *par);
        Double_t ExponentialDecay(Double_t *x, Double_t *par);

        Double_t ExponentialDecay_RiseSingle(Double_t *x, Double_t *par);
        Double_t ExponentialDecay_RiseSum(Double_t *x, Double_t *par);
        Double_t ExponentialDecay_ConvSingle(Double_t *x, Double_t *par);
        Double_t ExponentialDecay_ConvSum(Double_t *x, Double_t *par);

        int summation_order;
        Double_t detector_resolution;

    private:
 
};

#endif
