#include "TimeStampTest.h"

TimeStampTest::TimeStampTest():Tool(){}


bool TimeStampTest::Initialise(std::string configfile, DataModel &data)
{

    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

    return true;
}


bool TimeStampTest::Execute()
{
    if(m_data->TD.Stop==true)
    {
        std::cout <<"Stop has been called, skipping Timestamp calculator"<<std::endl;
        return true;
    }

    if(m_data->TD.EndOfRun==true)
    {
        std::cout <<"Run-End has been called, skipping Timestamp calculator and cleaning up"<<std::endl;
        Timestamp_Ch0 = m_data->TD.Timestamp_Map[0];
        Timestamp_Ch1 = m_data->TD.Timestamp_Map[1];
        Timestamp_Ch2 = m_data->TD.Timestamp_Map[2];
        Timestamp_Ch3 = m_data->TD.Timestamp_Map[3];
        Timestamp_Ch4 = m_data->TD.Timestamp_Map[4];
        Timestamp_Ch5 = m_data->TD.Timestamp_Map[5];
        Timestamp_Ch6 = m_data->TD.Timestamp_Map[6];
        Timestamp_Ch7 = m_data->TD.Timestamp_Map[7];
        Timestamp_Ch8 = m_data->TD.Timestamp_Map[8];
        Timestamp_Ch9 = m_data->TD.Timestamp_Map[9];
        Timestamp_Ch10 = m_data->TD.Timestamp_Map[10];
        Timestamp_Ch11 = m_data->TD.Timestamp_Map[11];
        Timestamp_Ch12 = m_data->TD.Timestamp_Map[12];
        Timestamp_Ch13 = m_data->TD.Timestamp_Map[13];
        Timestamp_Ch14 = m_data->TD.Timestamp_Map[14];
        Timestamp_Ch15 = m_data->TD.Timestamp_Map[15];

        m_data->TD.RootFile_Analysis->cd();
        m_data->TD.TTree_Analysis_Timestamp->Fill();
        m_data->TD.TTree_Analysis_Timestamp->Write();
        return true;
    }

    if(m_data->TD.NewRun==true)
    {
        InitRoot();
    }

    for(int i_channel: m_data->TD.ListOfChannels)
    {
        if(m_data->TD.PeakPositions[i_channel].size()<=0){continue;}

        if(i_channel>=0 && i_channel<=7)
        {
            trg_ch = 16;
            TriggerThresholdChannel = m_data->TD.ParsedMap_Header[trg_ch][5];
        }else if(i_channel>=8 && i_channel<=15)
        {
            trg_ch = 17;
            TriggerThresholdChannel = m_data->TD.ParsedMap_Header[trg_ch][5];
        }

        ChannelThresholdChannel = m_data->TD.ParsedMap_Header[i_channel][5];

        m_data->TD.Timestamp_Map[i_channel].push_back((ChannelThresholdChannel-TriggerThresholdChannel));
        std::cout << "Just wrote " << ChannelThresholdChannel << " for a trigger " << TriggerThresholdChannel << std::endl; 
    }

    return true;
}


bool TimeStampTest::Finalise()
{
    return true;
}


void TimeStampTest::InitRoot()
{
    m_data->TD.RootFile_Analysis->cd();
    m_data->TD.TTree_Analysis_Timestamp = new TTree("TTree_Analysis_Timestamp", "TTree_Analysis_Timestamp");
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch0", &Timestamp_Ch0);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch1", &Timestamp_Ch1);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch2", &Timestamp_Ch2);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch3", &Timestamp_Ch3);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch4", &Timestamp_Ch4);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch5", &Timestamp_Ch5);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch6", &Timestamp_Ch6);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch7", &Timestamp_Ch7);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch8", &Timestamp_Ch8);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch9", &Timestamp_Ch9);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch10", &Timestamp_Ch10);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch11", &Timestamp_Ch11);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch12", &Timestamp_Ch12);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch13", &Timestamp_Ch13);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch14", &Timestamp_Ch14);
    m_data->TD.TTree_Analysis_Timestamp->Branch("Timestamp_Ch15", &Timestamp_Ch15);  
}