#include <FitModels.h>

FitModels::FitModels()
{
    detector_resolution = 1;
    summation_order = 0;
}

FitModels::~FitModels()
{}

Double_t FitModels::Gaussian(Double_t *x, Double_t *par)
{
    Double_t amp = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];

    return (amp / (TMath::Sqrt(2.0 * TMath::Pi())*sigma)) * TMath::Exp(-TMath::Power((x[0]-mean),2)/(2*TMath::Power(sigma,2)));
}

Double_t FitModels::ExponentialDecay(Double_t *x, Double_t *par)
{
    Double_t n = par[0];
    Double_t tau = par[1];
    Double_t t0 = par[2];

    return (n/tau) * TMath::Exp(-(x[0]-t0)/tau);
}

Double_t FitModels::ExponentialDecay_RiseSingle(Double_t *x, Double_t *par)
{
    Double_t x0 = par[0];
    Double_t tau_R = par[1];
    Double_t n = par[2];
    Double_t tau = par[3];

    // Add if gauss convolution
    // 1/(TMath::Sqrt(2*TMath::Pi())*detector_resolution) * TMath::Exp(-TMath::Power(x[0]-x0,2) / (2*TMath::Power(detector_resolution,2))) *

    return n * (  (   TMath::Exp(-(x[0]-x0)/tau) - TMath::Exp(-(x[0]-x0)/tau_R) ) /   (tau - tau_R)   );
}

Double_t FitModels::ExponentialDecay_RiseSum(Double_t *x, Double_t *par)
{
    Double_t x0, tau_R, n_0, tau_0, n_1, tau_1, n_2, tau_2, n_3, tau_3;   

    if(summation_order>=1)
    {
        x0 = par[0];
        tau_R = par[1];
        n_0 = par[2];
        tau_0 = par[3];
    }
    if(summation_order>=2)
    {
        n_1 = par[4];
        tau_1 = par[5];
    }
    if(summation_order>=3)
    {
        n_2 = par[6];
        tau_2 = par[7];
    }
    if(summation_order>=4)
    {
        n_3 = par[8];
        tau_3 = par[9];
    }

    Double_t ReturnValue = n_0 * (  (   TMath::Exp(-(x[0]-x0)/tau_0) - TMath::Exp(-(x[0]-x0)/tau_R) ) /   (tau_0 - tau_R)   );
    if(summation_order>=2)
    {
        ReturnValue += n_1 * (  (   TMath::Exp(-(x[0]-x0)/tau_1) - TMath::Exp(-(x[0]-x0)/tau_R) ) /   (tau_1 - tau_R)   );
    }
    if(summation_order>=3)
    {
        ReturnValue += n_2 * (  (   TMath::Exp(-(x[0]-x0)/tau_2) - TMath::Exp(-(x[0]-x0)/tau_R) ) /   (tau_2 - tau_R)   );
    }
    if(summation_order>=4)
    {
        ReturnValue += n_3 * (  (   TMath::Exp(-(x[0]-x0)/tau_3) - TMath::Exp(-(x[0]-x0)/tau_R) ) /   (tau_3 - tau_R)   );
    }

    // Add if gauss convolution
    // 1/(TMath::Sqrt(2*TMath::Pi())*detector_resolution) * TMath::Exp(-TMath::Power(x[0]-x0,2) / (2*TMath::Power(detector_resolution,2))) *
    return ReturnValue;
}

Double_t FitModels::ExponentialDecay_ConvSingle(Double_t *x, Double_t *par)
{
    Double_t n = par[0];
    Double_t tau = par[1];
    Double_t x0 = par[2];

    return TMath::Exp(-TMath::Power(x[0]-x0,2) / (2*TMath::Power(detector_resolution,2))) * (n/tau) * TMath::Exp(-(x[0]-x0)/tau);
}

Double_t FitModels::ExponentialDecay_ConvSum(Double_t *x, Double_t *par)
{
    Double_t x0, n_0, tau_0, n_1, tau_1, n_2, tau_2, n_3, tau_3;   

    if(summation_order>=1)
    {
        x0 = par[0];
        n_0 = par[1];
        tau_0 = par[2];
    }
    if(summation_order>=2)
    {
        n_1 = par[3];
        tau_1 = par[4];
    }
    if(summation_order>=3)
    {
        n_2 = par[5];
        tau_2 = par[6];
    }
    if(summation_order>=4)
    {
        n_3 = par[7];
        tau_3 = par[8];
    }

    Double_t ReturnValue = (n_0/tau_0) * TMath::Exp(-(x[0] - x0)/tau_0);
    if(summation_order>=2)
    {
        ReturnValue += (n_1/tau_1) * TMath::Exp(-(x[0] - x0)/tau_1);
    }
    if(summation_order>=3)
    {
        ReturnValue += (n_2/tau_2) * TMath::Exp(-(x[0] - x0)/tau_2);
    }
    if(summation_order>=4)
    {
        ReturnValue += (n_3/tau_3) * TMath::Exp(-(x[0] - x0)/tau_3);
    }


    return TMath::Exp(-TMath::Power(x[0]-x0,2) / (2*TMath::Power(detector_resolution,2))) * ReturnValue;
}