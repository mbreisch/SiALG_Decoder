#include "LoadBinaryFile.h"

LoadBinaryFile::LoadBinaryFile():Tool(){}


bool LoadBinaryFile::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
    if(!m_variables.Get("Runtime",m_data->TD.Runtime)) m_data->TD.Runtime=-1;

    //Get Path input and filename
    if(!m_variables.Get("Path_In",m_data->TD.Path_In)) m_data->TD.Path_In="./";
    if(!m_variables.Get("File_Prefix",File_Prefix)) File_Prefix="wave_";
    if(!m_variables.Get("File_Sufix",File_Sufix)) File_Sufix=".dat";
    if(!m_variables.Get("Trigger_0",Trigger_0)) Trigger_0="TR_0_0.dat";
    if(!m_variables.Get("Trigger_1",Trigger_1)) Trigger_1="TR_0_1.dat";

    if(!m_variables.Get("words_in_data",words_in_data)) words_in_data=0;
    if(!m_variables.Get("bytes_in_dataword",bytes_in_dataword)) bytes_in_dataword=0;

    if(!m_variables.Get("words_in_header",words_in_header)) words_in_header=0;
    if(!m_variables.Get("bytes_in_headerword",bytes_in_headerword)) bytes_in_headerword=0;

    std::cout<<"Loading data files"<<std::endl;
    for(int i_channel=0; i_channel<MAX_NUM_CHANNELS; i_channel++)
    {
        FullFilePath = m_data->TD.Path_In + File_Prefix + std::to_string(i_channel) + File_Sufix;
        ifstream tmpfile(FullFilePath, std::ios::binary);
        if(!tmpfile) 
        {
            std::cout << "Channel " << i_channel << " not found" << std::endl;
            continue;
        }
        m_data->TD.FileMap.insert(std::make_pair(i_channel, std::move(tmpfile)));
        tmpfile.close();
    }
    std::cout<<"Loading trigger 0 files"<<std::endl;
    FullFilePath = m_data->TD.Path_In + Trigger_0;
    ifstream tmpfile2(FullFilePath, std::ios::binary);
    if(!tmpfile2) 
    {
        std::cout << "Channel Trigger 0 not found" << std::endl;
    }
    m_data->TD.FileMap.insert(std::make_pair(16, std::move(tmpfile2)));
    tmpfile2.close();
    std::cout<<"Loading trigger 1 files"<<std::endl;
    FullFilePath = m_data->TD.Path_In + Trigger_1;
    ifstream tmpfile3(FullFilePath, std::ios::binary);
    if(!tmpfile3) 
    {
        std::cout << "Channel Trigger 1 not found" << std::endl;
    }
    m_data->TD.FileMap.insert(std::make_pair(17, std::move(tmpfile3)));
    tmpfile3.close();

    std::cout << "Got " << m_data->TD.FileMap.size() << " to read out" << std::endl;

    m_data->TD.EventCounter = 0;
    
    return true;
}


bool LoadBinaryFile::Execute()
{
    m_data->TD.ParsedMap_Data.clear();
    m_data->TD.ParsedMap_Header.clear();
    if(m_data->TD.EventCounter>=m_data->TD.Runtime && m_data->TD.Runtime!=-1)
    {
        std::cout << "Eventcout = " << m_data->TD.EventCounter << " exceded the runtime of " << m_data->TD.Runtime << std::endl;
        m_data->TD.Stop = true;
        return true;
    }

    std::cout << "Starting with event " << m_data->TD.EventCounter << std::endl;  

    //Data
    for(int i_channel=0; i_channel<MAX_NUM_CHANNELS; i_channel++)
    {
        if(!m_data->TD.FileMap[i_channel].eof())
        {
            tmp_H = LoadHeader(i_channel);
            tmp_D = LoadData(i_channel);
            if(!tmp_D.empty() && !tmp_H.empty())
            {
            m_data->TD.ParsedMap_Header.insert(std::pair<int,vector<unsigned int>>(i_channel,tmp_H));
            m_data->TD.ParsedMap_Data.insert(std::pair<int,vector<float>>(i_channel,tmp_D));
            }
        }else
        {
            std::cout << "Reached EOF for channel " << i_channel << std::endl;
            m_data->TD.Stop = true;
        }
    }

    //Trigger 0
    int o_channel = 16;
    if(!m_data->TD.FileMap[o_channel].eof())
    {
        tmp_H = LoadHeader(o_channel);
        tmp_D = LoadData(o_channel);
        if(!tmp_D.empty() && !tmp_H.empty())
        {
            m_data->TD.ParsedMap_Header.insert(std::pair<int,vector<unsigned int>>(o_channel,tmp_H));
            m_data->TD.ParsedMap_Data.insert(std::pair<int,vector<float>>(o_channel,tmp_D));
        }
    }else
    {
        std::cout << "Reached EOF for channel " << o_channel << std::endl;
        m_data->TD.Stop = true;
    }

    //Trigger 1
    o_channel = 17;
    if(!m_data->TD.FileMap[o_channel].eof())
    {
        tmp_H = LoadHeader(o_channel);
        tmp_D = LoadData(o_channel);
        if(!tmp_D.empty() && !tmp_H.empty())
        {
            m_data->TD.ParsedMap_Header.insert(std::pair<int,vector<unsigned int>>(o_channel,tmp_H));
            m_data->TD.ParsedMap_Data.insert(std::pair<int,vector<float>>(o_channel,tmp_D));
        }
    }else
    {
        std::cout << "Reached EOF for channel " << o_channel << std::endl;
        m_data->TD.Stop = true;
    }

    if(!m_data->TD.Stop)
    {
        m_data->TD.EventCounter++;
    }

    return true;
}


bool LoadBinaryFile::Finalise()
{
    for(int i_channel=0; i_channel<MAX_NUM_CHANNELS; i_channel++)
    {
        m_data->TD.FileMap[i_channel].close();
    }
    m_data->TD.FileMap[16].close();
    m_data->TD.FileMap[17].close();
    m_data->TD.FileMap.clear();

    return true;
}

vector<unsigned int> LoadBinaryFile::LoadHeader(int i_channel)
{
    int buffersize = words_in_header * bytes_in_headerword;
    char buffer[buffersize];

    if (!m_data->TD.FileMap[i_channel].read(buffer, buffersize)) 
    {
        std::cerr << "Error occurred while reading the file for channel " << i_channel << std::endl;
        return {};
    }

    int numFloats = buffersize / sizeof(unsigned int);
    std::vector<unsigned int> Header(numFloats);
    std::memcpy(Header.data(), buffer, buffersize);

    if(m_verbose>2)
    {
        for(unsigned int k: Header)
        {
            std::cout<<"H :"<<k<<std::endl;
        }
    }
    return Header;
}

vector<float> LoadBinaryFile::LoadData(int i_channel)
{
    int buffersize = words_in_data * bytes_in_dataword;
    char buffer[buffersize];

    if(!m_data->TD.FileMap[i_channel].read(buffer, buffersize)) 
    {
        std::cerr << "Error occurred while reading the file." << std::endl;
        return {};
    }

    int numFloats = buffersize / sizeof(float);
    std::vector<float> Data(numFloats);
    std::memcpy(Data.data(), buffer, buffersize);

    if(m_verbose>2)
    {
        for(float k: Data)
        {
            std::cout<<"D :"<<k<<std::endl;
        }
    }

    return Data;   
}

