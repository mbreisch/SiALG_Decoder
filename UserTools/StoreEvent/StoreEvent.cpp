#include "StoreEvent.h"

StoreEvent::StoreEvent():Tool(){}


bool StoreEvent::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

    if(!m_variables.Get("Path_Out",m_data->TD.Path_Out)) m_data->TD.Path_Out="./";

    //Generate the rootfile
    std::string savelocation = m_data->TD.Path_Out+ "RawData.root";
    m_data->TD.RootFile_Event = new TFile(savelocation.c_str(),"RECREATE");
    m_data->TD.TTree_Event = new TTree("Event", "Event");
    m_data->TD.TTree_Event->Branch("EventID", &EventID, "EventID/I");

    m_data->TD.TTree_Event->Branch("Header_ch0", &Header_ch0);
    m_data->TD.TTree_Event->Branch("Data_ch0", &Data_ch0);

    m_data->TD.TTree_Event->Branch("Header_ch1", &Header_ch1);
    m_data->TD.TTree_Event->Branch("Data_ch1", &Data_ch1);

    m_data->TD.TTree_Event->Branch("Header_ch2", &Header_ch2);
    m_data->TD.TTree_Event->Branch("Data_ch2", &Data_ch2);
    
    m_data->TD.TTree_Event->Branch("Header_ch3", &Header_ch3);
    m_data->TD.TTree_Event->Branch("Data_ch3", &Data_ch3);

    m_data->TD.TTree_Event->Branch("Header_ch4", &Header_ch4);
    m_data->TD.TTree_Event->Branch("Data_ch4", &Data_ch4);

    m_data->TD.TTree_Event->Branch("Header_ch5", &Header_ch5);
    m_data->TD.TTree_Event->Branch("Data_ch5", &Data_ch5);

    m_data->TD.TTree_Event->Branch("Header_ch5", &Header_ch5);
    m_data->TD.TTree_Event->Branch("Data_ch5", &Data_ch5);
    
    m_data->TD.TTree_Event->Branch("Header_ch6", &Header_ch6);
    m_data->TD.TTree_Event->Branch("Data_ch6", &Data_ch6);
    
    m_data->TD.TTree_Event->Branch("Header_ch7", &Header_ch7);
    m_data->TD.TTree_Event->Branch("Data_ch7", &Data_ch7);

    m_data->TD.TTree_Event->Branch("Header_ch8", &Header_ch8);
    m_data->TD.TTree_Event->Branch("Data_ch8", &Data_ch8);

    m_data->TD.TTree_Event->Branch("Header_ch9", &Header_ch9);
    m_data->TD.TTree_Event->Branch("Data_ch9", &Data_ch9);
    
    m_data->TD.TTree_Event->Branch("Header_ch10", &Header_ch10);
    m_data->TD.TTree_Event->Branch("Data_ch10", &Data_ch10);

    m_data->TD.TTree_Event->Branch("Header_ch11", &Header_ch11);
    m_data->TD.TTree_Event->Branch("Data_ch11", &Data_ch11);

    m_data->TD.TTree_Event->Branch("Header_ch12", &Header_ch12);
    m_data->TD.TTree_Event->Branch("Data_ch12", &Data_ch12);

    m_data->TD.TTree_Event->Branch("Header_ch13", &Header_ch13);
    m_data->TD.TTree_Event->Branch("Data_ch13", &Data_ch13);
    
    m_data->TD.TTree_Event->Branch("Header_ch14", &Header_ch14);
    m_data->TD.TTree_Event->Branch("Data_ch14", &Data_ch14);

    m_data->TD.TTree_Event->Branch("Header_ch15", &Header_ch15);
    m_data->TD.TTree_Event->Branch("Data_ch15", &Data_ch15);

    m_data->TD.TTree_Event->Branch("Trigger_header_g0", &Trigger_header_g0);
    m_data->TD.TTree_Event->Branch("Trigger_data_g0", &Trigger_data_g0);

    m_data->TD.TTree_Event->Branch("Trigger_header_g1", &Trigger_header_g1);
    m_data->TD.TTree_Event->Branch("Trigger_data_g1", &Trigger_data_g1);

    return true;
}


bool StoreEvent::Execute()
{
    if(m_data->TD.Stop==true)
    {
        m_data->vars.Set("StopLoop",1);
        return true;
    }

    if(m_verbose>1){std::cout << "Got size " << m_data->TD.ParsedMap_Data.size() << std::endl;}
    //printMap(m_data->TD.ParsedMap);

    if(m_verbose>1){std::cout<<"Copying event ID..."<<std::endl;}
    EventID = m_data->TD.EventCounter;

    if(m_verbose>1){std::cout<<"Copying headers..."<<std::endl;}
    Header_ch0 = m_data->TD.ParsedMap_Header[0];
    Header_ch1 = m_data->TD.ParsedMap_Header[1];
    Header_ch2 = m_data->TD.ParsedMap_Header[2];
    Header_ch3 = m_data->TD.ParsedMap_Header[3];
    Header_ch4 = m_data->TD.ParsedMap_Header[4];
    Header_ch5 = m_data->TD.ParsedMap_Header[5];
    Header_ch6 = m_data->TD.ParsedMap_Header[6];
    Header_ch7 = m_data->TD.ParsedMap_Header[7];
    Header_ch8 = m_data->TD.ParsedMap_Header[8];
    Header_ch9 = m_data->TD.ParsedMap_Header[9];
    Header_ch10 = m_data->TD.ParsedMap_Header[10];
    Header_ch11 = m_data->TD.ParsedMap_Header[11];
    Header_ch12 = m_data->TD.ParsedMap_Header[12];
    Header_ch13 = m_data->TD.ParsedMap_Header[13];
    Header_ch14 = m_data->TD.ParsedMap_Header[14];
    Header_ch15 = m_data->TD.ParsedMap_Header[15];

    if(m_verbose>1){std::cout<<"Copying data..."<<std::endl;}
    Data_ch0 = m_data->TD.ParsedMap_Data[0];
    Data_ch1 = m_data->TD.ParsedMap_Data[1];
    Data_ch2 = m_data->TD.ParsedMap_Data[2];
    Data_ch3 = m_data->TD.ParsedMap_Data[3];
    Data_ch4 = m_data->TD.ParsedMap_Data[4];
    Data_ch5 = m_data->TD.ParsedMap_Data[5];
    Data_ch6 = m_data->TD.ParsedMap_Data[6];
    Data_ch7 = m_data->TD.ParsedMap_Data[7];
    Data_ch8 = m_data->TD.ParsedMap_Data[8];
    Data_ch9 = m_data->TD.ParsedMap_Data[9];
    Data_ch10 = m_data->TD.ParsedMap_Data[10];
    Data_ch11 = m_data->TD.ParsedMap_Data[11];
    Data_ch12 = m_data->TD.ParsedMap_Data[12];
    Data_ch13 = m_data->TD.ParsedMap_Data[13];
    Data_ch14 = m_data->TD.ParsedMap_Data[14];
    Data_ch15 = m_data->TD.ParsedMap_Data[15];

    if(m_verbose>1){std::cout<<"Copying trigger..."<<std::endl;}
    Trigger_header_g0 = m_data->TD.ParsedMap_Header[16];
    Trigger_header_g1 = m_data->TD.ParsedMap_Header[17];
    Trigger_data_g0 = m_data->TD.ParsedMap_Data[16];
    Trigger_data_g1 = m_data->TD.ParsedMap_Data[17];

    m_data->TD.RootFile_Event->cd();
    m_data->TD.TTree_Event->Fill();
    m_data->TD.TTree_Event->Write();
    m_data->TD.TTree_Event->Reset();

    return true;
}


bool StoreEvent::Finalise()
{
    m_data->TD.RootFile_Event->cd();
    m_data->TD.RootFile_Event->Close();
    delete m_data->TD.RootFile_Event;
    m_data->TD.RootFile_Event = 0;
    return true;
}


void StoreEvent::printMap(const std::map<int, std::vector<std::vector<unsigned int>>>& myMap)
{
    for (const auto& entry : myMap)
    {
        std::cout << "Key: " << entry.first << std::endl;
        std::cout << "Value:" << std::endl;

        const auto& vectorVector = entry.second;
        for (const auto& innerVector : vectorVector)
        {
            for (const auto& value : innerVector)
            {
                std::cout << value << " ";
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;
    }
}