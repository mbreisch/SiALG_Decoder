#include "EndChain.h"

EndChain::EndChain():Tool(){}


bool EndChain::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

    if(!m_variables.Get("Sum_Charge",Sum_Charge)) Sum_Charge=0;
    if(!m_variables.Get("Sum_TTS",Sum_TTS)) Sum_TTS=0;
    if(!m_variables.Get("Sum_PHD",Sum_PHD)) Sum_PHD=0;


    return true;
}


bool EndChain::Execute()
{   
    if(m_data->TD.Stop==true)
    {
        std::cout << "Stopping loop now!" << std::endl;
        m_data->vars.Set("StopLoop",1);
        return true;
    }

    if(m_data->TD.EndOfRun==true)
    {
        std::cout << "Finalising Event Rootfile: ... ";
        m_data->TD.RootFile_Event->cd();
        m_data->TD.RootFile_Event->Close();
        // delete m_data->TD.RootFile_Event;
        // m_data->TD.RootFile_Event = 0;
        std::cout << "Done!" << std::endl;

        std::cout << "Finalising Analysis Rootfile: ... ";
        m_data->TD.RootFile_Analysis->cd();
        m_data->TD.RootFile_Analysis->Close();
        // delete m_data->TD.RootFile_Analysis;
        // m_data->TD.RootFile_Analysis = 0;
        std::cout << "Done!" << std::endl;

        //Cleanup
        for(int i_channel: m_data->TD.ListOfChannels)
        {
            std::cout << "Closing file for channel " << i_channel << std::endl;
            if(m_data->TD.FileMap[i_channel].is_open())
            {
                m_data->TD.FileMap[i_channel].close();
            }
        }
        std::cout<<"-------------------"<<std::endl;
        m_data->TD.ListOfChannels.clear();
        m_data->TD.FileMap.clear();
        if(Sum_PHD==0){m_data->TD.PeakMinima.clear();}
        if(Sum_TTS==0){m_data->TD.TTS_Map.clear();}
        if(Sum_Charge==0){m_data->TD.Charge_Map.clear();}

        m_data->TD.NewRun = true;
        m_data->TD.EndOfRun = false;
        m_data->TD.EventCounter = 0;
        return true;
    }

    if(m_data->TD.NewRun == true)
    {
        m_data->TD.NewRun = false;
        m_data->TD.EndOfRun = false;
    }

    return true;
}


bool EndChain::Finalise()
{
    return true;
}
