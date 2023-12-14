#include "GetTTS.h"

GetTTS::GetTTS():Tool(){}


bool GetTTS::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
    if(!m_variables.Get("ROI_low",ROI_low)) ROI_low=0;
    if(!m_variables.Get("ROI_high",ROI_high)) ROI_high=1023;
    if(!m_variables.Get("threshold_multiplier",threshold_multiplier)) threshold_multiplier=0.5;      
    if(!m_variables.Get("MultiplicityCut",MultiplicityCut)) MultiplicityCut=1;    

    return true;
}


bool GetTTS::Execute()
{
    if(m_data->TD.Stop==true)
    {
        std::cout <<"Stop has been called, skipping TTS calculator"<<std::endl;
        return true;
    }

    if(m_data->TD.EndOfRun==true)
    {
        std::cout <<"Run-End has been called, skipping TTS calculator and cleaning up"<<std::endl;
        TTS_Ch0 = m_data->TD.TTS_Map[0];
        TTS_Ch1 = m_data->TD.TTS_Map[1];
        TTS_Ch2 = m_data->TD.TTS_Map[2];
        TTS_Ch3 = m_data->TD.TTS_Map[3];
        TTS_Ch4 = m_data->TD.TTS_Map[4];
        TTS_Ch5 = m_data->TD.TTS_Map[5];
        TTS_Ch6 = m_data->TD.TTS_Map[6];
        TTS_Ch7 = m_data->TD.TTS_Map[7];
        TTS_Ch8 = m_data->TD.TTS_Map[8];
        TTS_Ch9 = m_data->TD.TTS_Map[9];
        TTS_Ch10 = m_data->TD.TTS_Map[10];
        TTS_Ch11 = m_data->TD.TTS_Map[11];
        TTS_Ch12 = m_data->TD.TTS_Map[12];
        TTS_Ch13 = m_data->TD.TTS_Map[13];
        TTS_Ch14 = m_data->TD.TTS_Map[14];
        TTS_Ch15 = m_data->TD.TTS_Map[15];

        m_data->TD.RootFile_Analysis->cd();
        m_data->TD.TTree_Analysis_TTS->Fill();
        m_data->TD.TTree_Analysis_TTS->Write();
        return true;
    }

    if(m_data->TD.NewRun==true)
    {
        InitRoot();
    }

    std::map<int,vector<float>> Single_Event_Map;

    for(int i_channel: m_data->TD.ListOfChannels)
    {
        if(m_data->TD.PeakPositions[i_channel].size()<=0){continue;}
        if(i_channel==16 || i_channel==17){continue;}

        if(i_channel>=0 && i_channel<=7)
        {
            trg_ch = 16;
            TriggerThresholdChannel = GetTriggerThreshold(trg_ch);
        }else if(i_channel>=8 && i_channel<=15)
        {
            trg_ch = 17;
            TriggerThresholdChannel = GetTriggerThreshold(trg_ch);
        }

        for(int peak: m_data->TD.PeakPositions[i_channel])
        {
            if(peak>=ROI_low && peak<=ROI_high)
            {
                vector<float> data = m_data->TD.ParsedMap_Data[i_channel];
                int position = peak;
                float amplitude = TMath::MinElement(data.size(),data.data());

                ChannelThresholdChannel = GetChannelThreshold(data,position,amplitude);
                m_data->TD.TTS_Map[i_channel].push_back((ChannelThresholdChannel-TriggerThresholdChannel)/5.12);
                Single_Event_Map[i_channel].push_back((ChannelThresholdChannel-TriggerThresholdChannel)/5.12);
            }
        }
    }

    if(Single_Event_Map.size()>=MultiplicityCut)
    {
        for(const auto& pair : Single_Event_Map) 
        {
            int key = pair.first;
            const auto& vector2 = pair.second;

            // cout<<"Channel "<<key<<" has "<<vector2.size()<<" with multipicity "<<MultiplicityCut<<endl;

            auto it = m_data->TD.TTS_Map_cut.find(key);
            if(it != m_data->TD.TTS_Map_cut.end()) 
            {
                // Key exists in TTS_Map_cut, extend the vector
                it->second.insert(it->second.end(), vector2.begin(), vector2.end());
            }else 
            {
                // Key does not exist in TTS_Map_cut, add a new entry
                m_data->TD.TTS_Map_cut[key] = vector2;
            }
        }
    }

    return true;
}


bool GetTTS::Finalise()
{
    return true;
}


Double_t GaussianWA(Double_t *x, Double_t *par)
{
    Double_t amp = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];
    return (amp) * TMath::Exp(-0.5 * TMath::Power((x[0] - mean) / sigma, 2));
}

Double_t sigmoid(Double_t *x, Double_t *par) 
{
    Double_t t = x[0];
    Double_t x0 = par[0];
    Double_t s = par[1];
    Double_t a = par[2];
    Double_t y = 1.0 / (1.0 + TMath::Exp((x0 - t) / s));
    return a * y;
}

float GetTTS::GetTriggerThreshold(int trg_ch)
{
    // TApplication app("myApp", nullptr, nullptr);
    float TriggerThresholdChannel = -1;

    vector<float> data = m_data->TD.ParsedMap_Data[trg_ch];
    TGraph* graph = new TGraph();
    for(int i=0; i<data.size(); i++)
    {
        graph->SetPoint(i, i, data[i]);     
    }

    float min_val = TMath::MinElement(data.size(),data.data());
    float min_pos = 0;
    for(int i=0; i<data.size(); i++)
    {
        if(min_val==data.at(i))
        {
            min_pos = i;
            break;
        }
    } 

    // Create a TF1 object using the fitting function and set initial parameter values
    TF1* fitFunc = new TF1("fitFunc", sigmoid, min_pos-50, min_pos+5, 3);
    fitFunc->SetParameters(min_pos,3,min_val); // Set initial parameter values

    // Fit the graph with the TF1 function
    graph->Fit(fitFunc, "RQ");

    // Access the fit parameters and their uncertainties
    Double_t r_x0 = fitFunc->GetParameter(0);
    Double_t r_s0 = fitFunc->GetParameter(1);
    Double_t r_a0 = fitFunc->GetParameter(2);

    threshold_channel = threshold_multiplier*min_val;

    // Double_t amplitudeError = fitFunc->GetParError(0);
    // Double_t meanError = fitFunc->GetParError(1);
    // Double_t stddevError = fitFunc->GetParError(2);

    TriggerThresholdChannel = r_x0 - r_s0 * TMath::Log(r_a0/threshold_channel-1);


    // Visualize the fit result
    // TCanvas* canvas = new TCanvas("canvas", "Fit Result");
    // graph->Draw("AP"); // Draw the original graph
    // fitFunc->Draw("same"); // Draw the fitted function on top
    // canvas->Update();
    // std::string name = m_data->TD.Path_Out+"Event_"+to_string(m_data->TD.EventCounter)+".png";
    // canvas->SaveAs(name.c_str());  // Save the plot to a file
    // canvas->Draw();
    // app.Run();

    delete graph;
    delete fitFunc;

    return TriggerThresholdChannel;
}


float GetTTS::GetChannelThreshold(vector<float> data, int position,float amplitude)
{
    // TApplication app("myApp", nullptr, nullptr);
    float ChannelThresholdChannel = -1;

    TGraph* graph = new TGraph();
    for(int i=0; i<data.size(); i++)
    {
        graph->SetPoint(i, i, data[i]);     
    }

    // Create a TF1 object using the fitting function and set initial parameter values
    TF1* fitFunc = new TF1("fitFunc", GaussianWA, ROI_low, ROI_high, 3);
    fitFunc->SetParameters(amplitude, position, 5.0); // Set initial parameter values

    // Fit the graph with the TF1 function
    graph->Fit(fitFunc, "RQ");

    // Access the fit parameters and their uncertainties
    Double_t r_amplitude = fitFunc->GetParameter(0);
    Double_t r_mean = fitFunc->GetParameter(1);
    Double_t r_sigma = fitFunc->GetParameter(2);

    threshold_channel = threshold_multiplier*r_amplitude;

    // Double_t amplitudeError = fitFunc->GetParError(0);
    // Double_t meanError = fitFunc->GetParError(1);
    // Double_t stddevError = fitFunc->GetParError(2);

    float t1 = r_mean + TMath::Sqrt( - 2 * TMath::Power(r_sigma,2) * TMath::Log(threshold_channel/r_amplitude));
    float t2 = r_mean - TMath::Sqrt( - 2 * TMath::Power(r_sigma,2) * TMath::Log(threshold_channel/r_amplitude));

    if(t1>t2)
    {
        ChannelThresholdChannel = t2;
    }else
    {
        ChannelThresholdChannel = t1;
    }

    // Visualize the fit result
    // TCanvas* canvas = new TCanvas("canvas", "Fit Result");
    // graph->Draw("AP"); // Draw the original graph
    // fitFunc->Draw("same"); // Draw the fitted function on top
    // canvas->Draw();
    // app.Run();

    delete graph;
    delete fitFunc;

    return ChannelThresholdChannel;
}


void GetTTS::InitRoot()
{
    m_data->TD.RootFile_Analysis->cd();
    m_data->TD.TTree_Analysis_TTS = new TTree("TTree_Analysis_TTS", "TTree_Analysis_TTS");
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch0", &TTS_Ch0);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch1", &TTS_Ch1);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch2", &TTS_Ch2);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch3", &TTS_Ch3);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch4", &TTS_Ch4);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch5", &TTS_Ch5);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch6", &TTS_Ch6);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch7", &TTS_Ch7);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch8", &TTS_Ch8);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch9", &TTS_Ch9);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch10", &TTS_Ch10);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch11", &TTS_Ch11);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch12", &TTS_Ch12);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch13", &TTS_Ch13);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch14", &TTS_Ch14);
    m_data->TD.TTree_Analysis_TTS->Branch("TTS_Ch15", &TTS_Ch15);  
}