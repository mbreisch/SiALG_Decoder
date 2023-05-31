#include "GetEvent.h"

GetEvent::GetEvent():Tool(){}


bool GetEvent::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;


    return true;
}


bool GetEvent::Execute()
{
    return true;
}


bool GetEvent::Finalise()
{
    return true;
}
