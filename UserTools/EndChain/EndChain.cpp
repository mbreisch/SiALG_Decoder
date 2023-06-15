#include "EndChain.h"

EndChain::EndChain():Tool(){}


bool EndChain::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;


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

    return true;
}


bool EndChain::Finalise()
{
    std::cout << "Finalising Event Rootfile: ... ";
    m_data->TD.RootFile_Event->cd();
    m_data->TD.RootFile_Event->Close();
    delete m_data->TD.RootFile_Event;
    m_data->TD.RootFile_Event = 0;
    std::cout << "Done!" << std::endl;

    std::cout << "Finalising Analysis Rootfile: ... ";
    m_data->TD.RootFile_Analysis->cd();
    m_data->TD.RootFile_Analysis->Close();
    delete m_data->TD.RootFile_Analysis;
    m_data->TD.RootFile_Analysis = 0;
    std::cout << "Done!" << std::endl;
    return true;
}
