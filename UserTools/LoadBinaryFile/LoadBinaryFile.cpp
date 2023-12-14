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

    m_data->TD.EventCounter = 0;
    TMP_InPath = "startline";

    RunList.open(m_data->TD.Path_In,ios_base::in); // Open the file
    if(!RunList.is_open()) 
    {
        std::cerr << "Error opening RunList." << std::endl;
        return 1;
    }

    return true;
}


bool LoadBinaryFile::Execute()
{
    if(m_data->TD.Stop==true)
    {
        std::cout <<"Stop should never be called here!!!"<<std::endl;
        return true;
    }

    m_data->TD.ParsedMap_Data.clear();
    m_data->TD.ParsedMap_Header.clear();
    if(m_data->TD.EventCounter>=m_data->TD.Runtime && m_data->TD.Runtime!=-1)
    {
        std::cout << "Eventcount=" << m_data->TD.EventCounter << " reached condition of " << m_data->TD.Runtime << " events, called stop." << std::endl;
        m_data->TD.EndOfRun = true;
        return true;
    }

    if(m_data->TD.NewRun==true)
    {
        if(!LoadNewRun())
        {
            return true;
        }
    }

    if(m_verbose>1){std::cout << "-------------------------" << std::endl;}
    if(m_verbose>1){std::cout << "Starting with event " << m_data->TD.EventCounter+1 << std::endl;}
    else{std::cout << "\rStarting with event " << m_data->TD.EventCounter+1 << "                   " << std::flush;}

    //Data
    for(int i_channel: m_data->TD.ListOfChannels)
    {
        if(!m_data->TD.FileMap[i_channel].eof())
        {
            tmp_H = LoadHeader(i_channel);
            tmp_D = LoadData(i_channel);
            if(!tmp_D.empty() && !tmp_H.empty())
            {
                m_data->TD.ParsedMap_Header.insert(std::pair<int,vector<unsigned int>>(i_channel,tmp_H));
                m_data->TD.ParsedMap_Data.insert(std::pair<int,vector<float>>(i_channel,tmp_D));
            }else
            {
                
                if(m_verbose>2)
                {
                    std::cout<<"-------------------"<<std::endl;
                    std::cout << "Reached EOF 1 for channel " << i_channel << std::endl;
                }
                m_data->TD.EndOfRun = true;
            }
        }else
        {
            std::cout << "Reached EOF 2 for channel " << i_channel << std::endl;
            m_data->TD.EndOfRun = true;
        }
    }

    if(!m_data->TD.EndOfRun)
    {
        m_data->TD.EventCounter++;
    }

    return true;
}


bool LoadBinaryFile::Finalise()
{
    RunList.close();
    return true;
}

vector<unsigned int> LoadBinaryFile::LoadHeader(int i_channel)
{
    int buffersize = words_in_header * bytes_in_headerword;
    char buffer[buffersize];

    m_data->TD.FileMap[i_channel].read(buffer, buffersize);
    if(m_data->TD.FileMap[i_channel].fail()) 
    {
        return {};
    }else if(m_data->TD.FileMap[i_channel].bad())
    {
        std::cerr << "Non-recoverable read error occurred while reading the file for channel " << i_channel << std::endl;
        return {};
    }else if(m_data->TD.FileMap[i_channel].eof())
    {
        if (m_data->TD.FileMap[i_channel].gcount() == 0) 
        {
            std::cerr << "Premature EOF reading the file for channel " << i_channel << std::endl;
            return {};
        }else 
        {
            std::cerr << "EOF reading the file for channel " << i_channel << " after " << m_data->TD.FileMap[i_channel].gcount() << " bytes." << std::endl;
            return {};
        }
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

    m_data->TD.FileMap[i_channel].read(buffer, buffersize);
    if(m_data->TD.FileMap[i_channel].fail()) 
    {
        return {};
    }else if(m_data->TD.FileMap[i_channel].bad())
    {
        std::cerr << "Non-recoverable read error occurred while reading the file for channel " << i_channel << std::endl;
        return {};
    }else if(m_data->TD.FileMap[i_channel].eof())
    {
        if (m_data->TD.FileMap[i_channel].gcount() == 0) 
        {
            std::cerr << "Premature EOF reading the file for channel " << i_channel << std::endl;
            return {};
        }else 
        {
            std::cerr << "EOF reading the file for channel " << i_channel << " after " << m_data->TD.FileMap[i_channel].gcount() << " bytes." << std::endl;
            return {};
        }
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


bool LoadBinaryFile::LoadNewRun()
{
    std::getline(RunList, line);

    while(true)
    {
        if(line == TMP_InPath || line.empty() == true) //TMP_InPath is "startline" as init value
        {
            m_data->TD.Stop = true;
            return false;
        }else if(line[0] == '#')
        {   
            std::getline(RunList, line); //if # is at the start of a line ignore it and get the next line
        }else
        {
            TMP_InPath = line;
            break;
        }
    }

    while(!TMP_InPath.empty() && (TMP_InPath.back() == '\n' || TMP_InPath.back() == '\r')) 
    {
        TMP_InPath.pop_back();
    }

    if(!m_variables.Get("Path_Out",m_data->TD.Path_Out)) m_data->TD.Path_Out=TMP_InPath;

    std::cout<<"Loading data files from "<<TMP_InPath<<std::endl;

    std::string logbook_name = m_data->TD.Path_In + "_logbook.txt";
    ofstream outfile_logbook(logbook_name.c_str(), ios_base::out | ios_base::app);
    outfile_logbook<<"Loading data files from "<<TMP_InPath<<std::endl;
    outfile_logbook.close();

    std::cout<<"-------------------"<<std::endl;
    for(int i_channel=0; i_channel<MAX_NUM_CHANNELS; i_channel++)
    {
        FullFilePath = TMP_InPath + File_Prefix + std::to_string(i_channel) + File_Sufix;
        ifstream tmpfile(FullFilePath, std::ios::binary);
        if(!tmpfile) 
        {
            std::cout << "\033[1;31m>> Channel " << i_channel << " not found\033[0m" << std::endl;
            continue;
        }
        std::cout << "\033[1;32m>> Channel " << i_channel << " found\033[0m" << std::endl;
        m_data->TD.ListOfChannels.push_back(i_channel);
        m_data->TD.FileMap.insert(std::make_pair(i_channel, std::move(tmpfile)));
        tmpfile.close();
    }

    FullFilePath = TMP_InPath + Trigger_0;
    ifstream tmpfile2(FullFilePath, std::ios::binary);
    if(!tmpfile2) 
    {
        std::cout << "\033[1;31m>> Channel Trigger 0 not found\033[0m" << std::endl;
    }else
    {
        std::cout << "\033[1;32m>> Channel Trigger 0 found\033[0m" << std::endl;
        m_data->TD.ListOfChannels.push_back(16);
        m_data->TD.FileMap.insert(std::make_pair(16, std::move(tmpfile2)));
        tmpfile2.close();
    }

    FullFilePath = TMP_InPath + Trigger_1;
    ifstream tmpfile3(FullFilePath, std::ios::binary);
    if(!tmpfile3) 
    {
        std::cout << "\033[1;31m>> Channel Trigger 1 not found\033[0m" << std::endl;
    }else
    {
        std::cout << "\033[1;32m>> Channel Trigger 1 found\033[0m" << std::endl;
        m_data->TD.ListOfChannels.push_back(17);
        m_data->TD.FileMap.insert(std::make_pair(17, std::move(tmpfile3)));
        tmpfile3.close();
    }
    std::cout<<"-------------------"<<std::endl;
    std::cout << "Got " << m_data->TD.ListOfChannels.size() << " channels to read out" << std::endl;

    std::string savelocation = m_data->TD.Path_Out+ "Analysis.root";
    m_data->TD.RootFile_Analysis = new TFile(savelocation.c_str(),"RECREATE");
    
    savelocation = m_data->TD.Path_Out+ "RawData.root";
    m_data->TD.RootFile_Event = new TFile(savelocation.c_str(),"RECREATE");
    
    return true;
}