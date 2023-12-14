#include "PulseIdentifier.h"

PulseIdentifier::PulseIdentifier():Tool(){}


bool PulseIdentifier::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
    if(!m_variables.Get("threshold",m_threshold)) m_threshold=2.0;
    if(!m_variables.Get("width",m_width)) m_width=5.0;
    if(!m_variables.Get("limit_to_first_hit",limit_to_first_hit)) limit_to_first_hit=0;
    

    return true;
}


bool PulseIdentifier::Execute()
{
    if(m_data->TD.Stop==true)
    {
        std::cout <<"Stop has been called, skipping Peakfinder"<<std::endl;
        return true;
    }

    if(m_data->TD.EndOfRun==true)
    {
        std::cout <<"Run-End has been called, skipping Peakfinder and cleaning up"<<std::endl;
        return true; 
    }

    if(m_threshold==0)
    {
        float old_threshold = m_threshold; 
        float voltage = GetCurrentVoltage();
        m_threshold = -5.0/3.25 * voltage + 457.0/13.0;
        if(old_threshold!=m_threshold)
        {
            std::cout<< "Old threshold was " << old_threshold << " and is now set to " << m_threshold <<std::endl;
        }
    }
    
    m_data->TD.PeakPositions.clear();
    vector<int> input;
    for(std::map<int,vector<float>>::iterator it=m_data->TD.ParsedMap_Data.begin(); it!=m_data->TD.ParsedMap_Data.end(); ++it)
    {
        if(it->first==16 || it->first==17)
        {
            input = FindWall(it->second);
            m_data->TD.PeakPositions.insert(std::pair<int,vector<int>>(it->first,input));
        }else
        {
            input = FindPulses(it->first,it->second);
            m_data->TD.PeakPositions.insert(std::pair<int,vector<int>>(it->first,input));
        }
    }
    return true;
}


bool PulseIdentifier::Finalise()
{
    return true;
}


float PulseIdentifier::GetCurrentVoltage()
{
    float voltage = 0.0;

    std::istringstream iss(m_data->TD.Path_Out.c_str());
    std::string token;
    std::vector<std::string> tokens;

    while (std::getline(iss, token, '/')) 
    {
        tokens.push_back(token);
    }

    // Get the last token
    if (!tokens.empty()) 
    {
        std::string lastToken = tokens.back();
        // Remove the last character
        if (!lastToken.empty()) 
        {
            lastToken.pop_back();
            try 
            {
                // Convert the modified token to a float
                size_t pos;
                float floatValue = std::stof(lastToken, &pos);

                // Check if the entire token is a valid float
                if (pos == lastToken.length()) 
                {
                    voltage = floatValue;
                } else {
                    std::cerr << "Invalid characters after the float in token: " << lastToken << std::endl;
                }
            }catch(const std::invalid_argument& e) 
            {
                // Handle the case where the conversion fails
                std::cerr << "Error converting to float for token '" << lastToken << "': " << e.what() << std::endl;
            }
        }else 
        {
            std::cerr << "Last token is empty." << std::endl;
        }
    }else 
    {
        std::cerr << "No tokens found." << std::endl;
    }

    return voltage;
}

std::vector<int> PulseIdentifier::FindPulses(int channel, std::vector<float> waveform) 
{
    float current_sample = 0.0;
    int waveform_length=0;
    bool pulsestarted=false;
    int low = 0;
    int high = 0;
    float minimum=0.0;
    int minimumpeak=-1;
    

    vector<int> PulseReturn = {};

    for(int i_sample=0;i_sample<waveform.size();i_sample++) 
    {
        current_sample = waveform.at(i_sample);

        if(current_sample<m_threshold) 
        {
            waveform_length++;
            if(!pulsestarted){low=i_sample;}
            pulsestarted=true;
        }else 
        {
            if(waveform_length<m_width) 
            {
                waveform_length=0; 
                pulsestarted=false; 
                low=0; 
                high=0;
            }else 
            {
                high=i_sample;
                for(int i=low; i<=high; i++) 
                {
                    float entry = waveform.at(i);
                    if(entry < minimum) 
                    {
                        minimum = entry;
                        minimumpeak = i;
                    }
                }
                PulseReturn.push_back(minimumpeak);
                if(m_verbose>1){std::cout<<"Pulse of "<<channel<<" started at "<<low<<" and ended at "<<high<<" , Position of peak was determined as "<<minimumpeak<<std::endl;}
                pulsestarted=false;
                waveform_length=0; 
                low=0; 
                high=0;
            }
        }
    }
    
    if(PulseReturn.size()>=limit_to_first_hit && limit_to_first_hit != -1)
    {
        vector<int> tmp_vec = PulseReturn;
        PulseReturn.clear();
        PulseReturn.insert(PulseReturn.begin(), tmp_vec.begin(), tmp_vec.begin() + limit_to_first_hit);
    }

    return PulseReturn;
}


std::vector<int> PulseIdentifier::FindWall(std::vector<float> waveform) 
{
    vector<int> PulseReturn;
    float threshold = -400.0;

    for(int i_sample=0;i_sample<waveform.size();i_sample++) 
    {
        if(waveform.at(i_sample)<-400)
        {
            PulseReturn.push_back(i_sample);
            break;
        }
    }

    return PulseReturn;
}