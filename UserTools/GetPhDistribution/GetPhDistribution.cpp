#include "GetPhDistribution.h"

GetPhDistribution::GetPhDistribution():Tool(){}


bool GetPhDistribution::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
    if(!m_variables.Get("ROI_low",ROI_low)) ROI_low=0;
    if(!m_variables.Get("ROI_high",ROI_high)) ROI_high=1023;


    m_data->TD.RootFile_Analysis->cd();
    m_data->TD.TTree_Analysis_PHD = new TTree("TTree_Analysis_PHD", "TTree_Analysis_PHD");
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch0", &PulseHeight_Ch0);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch1", &PulseHeight_Ch1);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch2", &PulseHeight_Ch2);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch3", &PulseHeight_Ch3);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch4", &PulseHeight_Ch4);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch5", &PulseHeight_Ch5);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch6", &PulseHeight_Ch6);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch7", &PulseHeight_Ch7);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch8", &PulseHeight_Ch8);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch9", &PulseHeight_Ch9);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch10", &PulseHeight_Ch10);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch11", &PulseHeight_Ch11);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch12", &PulseHeight_Ch12);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch13", &PulseHeight_Ch13);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch14", &PulseHeight_Ch14);
    m_data->TD.TTree_Analysis_PHD->Branch("PulseHeight_Ch15", &PulseHeight_Ch15);        

    return true;
}


bool GetPhDistribution::Execute()
{
    if(m_data->TD.Stop==true)
    {
        std::cout <<"Stop has been called, skipping Pulseheight calculator"<<std::endl;
        return true;
    }

    for(int i_channel: m_data->TD.ListOfChannels)
    {
        vector<float> InVector = GetSlicedDataFromROI(m_data->TD.ParsedMap_Data[i_channel],ROI_low,ROI_high);
        if(i_channel==16 || i_channel==17){continue;}
        if(m_data->TD.PeakPositions[i_channel].empty())
        {
            m_data->TD.PeakMinima[i_channel].push_back(TMath::MinElement(InVector.size(),InVector.data()));     
        }else
        {
            for(int pos: m_data->TD.PeakPositions[i_channel])
            {
                if(pos>=ROI_low && pos<=ROI_high)
                {
                    m_data->TD.PeakMinima[i_channel].push_back(m_data->TD.ParsedMap_Data[i_channel][pos]);        
                }
            }
        }
    }

    return true;
}


bool GetPhDistribution::Finalise()
{
    PulseHeight_Ch0 = m_data->TD.PeakMinima[0];
    PulseHeight_Ch1 = m_data->TD.PeakMinima[1];
    PulseHeight_Ch2 = m_data->TD.PeakMinima[2];
    PulseHeight_Ch3 = m_data->TD.PeakMinima[3];
    PulseHeight_Ch4 = m_data->TD.PeakMinima[4];
    PulseHeight_Ch5 = m_data->TD.PeakMinima[5];
    PulseHeight_Ch6 = m_data->TD.PeakMinima[6];
    PulseHeight_Ch7 = m_data->TD.PeakMinima[7];
    PulseHeight_Ch8 = m_data->TD.PeakMinima[8];
    PulseHeight_Ch9 = m_data->TD.PeakMinima[9];
    PulseHeight_Ch10 = m_data->TD.PeakMinima[10];
    PulseHeight_Ch11 = m_data->TD.PeakMinima[11];
    PulseHeight_Ch12 = m_data->TD.PeakMinima[12];
    PulseHeight_Ch13 = m_data->TD.PeakMinima[13];
    PulseHeight_Ch14 = m_data->TD.PeakMinima[14];
    PulseHeight_Ch15 = m_data->TD.PeakMinima[15];

    m_data->TD.RootFile_Analysis->cd();
    m_data->TD.TTree_Analysis_PHD->Fill();
    m_data->TD.TTree_Analysis_PHD->Write();

    return true;
}


std::vector<float> GetPhDistribution::GetSlicedDataFromROI(vector<float> data, int startIndex, int endIndex)
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
