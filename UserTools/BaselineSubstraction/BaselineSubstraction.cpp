#include "BaselineSubstraction.h"

BaselineSubstraction::BaselineSubstraction():Tool(){}


bool BaselineSubstraction::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
    if(!m_variables.Get("words_in_data",words_in_data)) words_in_data=0;

    return true;
}


bool BaselineSubstraction::Execute()
{
    for(std::map<int,vector<float>>::iterator it=m_data->TD.ParsedMap_Data.begin(); it!=m_data->TD.ParsedMap_Data.end(); ++it)
    {
        vector<float> in_vector;
        float mean = 0.0;
        if(it->first==16 || it->first==17)
        {
            in_vector = it->second;
            for(auto& entry : in_vector)
            {
                entry = 2*entry;
            }
            vector<float> SlicedVector = GetSlice(in_vector, 10, words_in_data/5);
            mean = TMath::Mean(SlicedVector.size(), SlicedVector.data());
        }else
        {
            in_vector = it->second;
            mean = FitGaussianAndGetMean(in_vector);
        }
        for(int i_entry=0; i_entry<in_vector.size(); i_entry++)
        {
            if(i_entry<10 || i_entry>words_in_data-10)
            {
                in_vector.at(i_entry) = 0;
            }else
            {
                in_vector.at(i_entry) -= mean;
            }
        }
        it->second = in_vector;
    }

    return true;
}


bool BaselineSubstraction::Finalise()
{
    return true;
}


Double_t Gaussian(Double_t *x, Double_t *par)
{
    Double_t mean = par[0];
    Double_t sigma = par[1];
    return (1.0 / (sigma * TMath::Sqrt(2.0 * TMath::Pi()))) * TMath::Exp(-0.5 * TMath::Power((x[0] - mean) / sigma, 2));
}


float BaselineSubstraction::FitGaussianAndGetMean(const std::vector<float>& data)
{
    // Create a histogram from the data
    Double_t min = *min_element(data.begin(), data.end());
    Double_t max = *max_element(data.begin(), data.end());
    TH1F *hist = new TH1F("hist", "Data Distribution", 50, min, max);

    for (const auto& value : data)
    {
        hist->Fill(value);
    }

    // Fit the Gaussian distribution
    TF1 *fitFunc = new TF1("fitFunc", Gaussian, min, max, 2);
    fitFunc->SetParameters(hist->GetMean(), hist->GetRMS());
    hist->Fit(fitFunc, "RNQ");

    // Extract the mean
    Double_t mean = fitFunc->GetParameter(0);

    delete hist;
    delete fitFunc;

    return (float)mean;
}


std::vector<float> BaselineSubstraction::GetSlice(vector<float> data, int startIndex, int endIndex)
{
    if (startIndex < 0 || endIndex > data.size() || startIndex >= endIndex)
    {
        // Invalid range, return an empty vector or handle the error as desired
        return vector<float>();
    }

    // Create a new vector containing the desired slice
    std::vector<float> slice(data.begin() + startIndex, data.begin() + endIndex);

    return slice;
}