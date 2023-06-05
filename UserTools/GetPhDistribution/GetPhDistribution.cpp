#include "GetPhDistribution.h"

GetPhDistribution::GetPhDistribution():Tool(){}


bool GetPhDistribution::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
    if(!m_variables.Get("Path_Out",m_data->TD.Path_Out)) m_data->TD.Path_Out="./";

    std::string savelocation = m_data->TD.Path_Out+ "Analysis.root";
    m_data->TD.RootFile_Analysis = new TFile(savelocation.c_str(),"RECREATE");
    m_data->TD.RootFile_Analysis->cd();
    m_data->TD.TTree_Analysis_PHD = new TTree("TTree_Analysis_PHD", "TTree_Analysis_PHD");
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch0", &PusleHeight_Ch0);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch1", &PusleHeight_Ch1);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch2", &PusleHeight_Ch2);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch3", &PusleHeight_Ch3);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch4", &PusleHeight_Ch4);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch5", &PusleHeight_Ch5);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch6", &PusleHeight_Ch6);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch7", &PusleHeight_Ch7);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch8", &PusleHeight_Ch8);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch9", &PusleHeight_Ch9);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch10", &PusleHeight_Ch10);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch11", &PusleHeight_Ch11);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch12", &PusleHeight_Ch12);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch13", &PusleHeight_Ch13);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch14", &PusleHeight_Ch14);
    m_data->TD.TTree_Analysis_PHD->Branch("PusleHeight_Ch15", &PusleHeight_Ch15);        

    for(int i_channel: m_data->TD.ListOfChannels)
    {
        if(i_channel==16 || i_channel==17){continue;}
        m_data->TD.PeakMinima.insert(std::pair<int,vector<float>>(i_channel,{}));
    }
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
        if(i_channel==16 || i_channel==17){continue;}
        if(m_data->TD.PeakPositions[i_channel].empty())
        {
            m_data->TD.PeakMinima[i_channel].push_back(TMath::Mean(m_data->TD.ParsedMap_Data[i_channel].size(),m_data->TD.ParsedMap_Data[i_channel].data()));     
        }else
        {
            for(int pos: m_data->TD.PeakPositions[i_channel])
            {
                m_data->TD.PeakMinima[i_channel].push_back(m_data->TD.ParsedMap_Data[i_channel][pos]);        
            }
        }
    }

    return true;
}


bool GetPhDistribution::Finalise()
{
    PusleHeight_Ch0 = m_data->TD.PeakMinima[0];
    PusleHeight_Ch1 = m_data->TD.PeakMinima[1];
    PusleHeight_Ch2 = m_data->TD.PeakMinima[2];
    PusleHeight_Ch3 = m_data->TD.PeakMinima[3];
    PusleHeight_Ch4 = m_data->TD.PeakMinima[4];
    PusleHeight_Ch5 = m_data->TD.PeakMinima[5];
    PusleHeight_Ch6 = m_data->TD.PeakMinima[6];
    PusleHeight_Ch7 = m_data->TD.PeakMinima[7];
    PusleHeight_Ch8 = m_data->TD.PeakMinima[8];
    PusleHeight_Ch9 = m_data->TD.PeakMinima[9];
    PusleHeight_Ch10 = m_data->TD.PeakMinima[10];
    PusleHeight_Ch11 = m_data->TD.PeakMinima[11];
    PusleHeight_Ch12 = m_data->TD.PeakMinima[12];
    PusleHeight_Ch13 = m_data->TD.PeakMinima[13];
    PusleHeight_Ch14 = m_data->TD.PeakMinima[14];
    PusleHeight_Ch15 = m_data->TD.PeakMinima[15];

    m_data->TD.RootFile_Analysis->cd();
    m_data->TD.TTree_Analysis_PHD->Fill();
    m_data->TD.TTree_Analysis_PHD->Write();

    m_data->TD.RootFile_Analysis->Close();
    delete m_data->TD.RootFile_Analysis;
    m_data->TD.RootFile_Analysis = 0;

    return true;
}
