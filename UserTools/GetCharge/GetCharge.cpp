#include "GetCharge.h"

GetCharge::GetCharge():Tool(){}


bool GetCharge::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
    if(!m_variables.Get("Visualisation",vis)) vis=0;
    if(!m_variables.Get("ROI_low",ROI_low)) ROI_low=0;
    if(!m_variables.Get("ROI_high",ROI_high)) ROI_high=1023;

    return true;
}


bool GetCharge::Execute()
{
    if(m_data->TD.Stop==true)
    {
        std::cout <<"Stop has been called, skipping Charge calculator"<<std::endl;
        return true;
    }

    if(m_data->TD.EndOfRun==true)
    {
        std::cout <<"Run-End has been called, skipping Charge calculator and cleaning up"<<std::endl;
        Charge_Ch0 = m_data->TD.Charge_Map[0];
        Charge_Ch1 = m_data->TD.Charge_Map[1];
        Charge_Ch2 = m_data->TD.Charge_Map[2];
        Charge_Ch3 = m_data->TD.Charge_Map[3];
        Charge_Ch4 = m_data->TD.Charge_Map[4];
        Charge_Ch5 = m_data->TD.Charge_Map[5];
        Charge_Ch6 = m_data->TD.Charge_Map[6];
        Charge_Ch7 = m_data->TD.Charge_Map[7];
        Charge_Ch8 = m_data->TD.Charge_Map[8];
        Charge_Ch9 = m_data->TD.Charge_Map[9];
        Charge_Ch10 = m_data->TD.Charge_Map[10];
        Charge_Ch11 = m_data->TD.Charge_Map[11];
        Charge_Ch12 = m_data->TD.Charge_Map[12];
        Charge_Ch13 = m_data->TD.Charge_Map[13];
        Charge_Ch14 = m_data->TD.Charge_Map[14];
        Charge_Ch15 = m_data->TD.Charge_Map[15];

        m_data->TD.RootFile_Analysis->cd();
        m_data->TD.TTree_Analysis_Charge->Fill();
        m_data->TD.TTree_Analysis_Charge->Write();
        return true;
    }

    if(m_data->TD.NewRun==true)
    {
        InitRoot();
    }

    for(int i_channel: m_data->TD.ListOfChannels)
    {
        vector<float> data = m_data->TD.ParsedMap_Data[i_channel];
        data = GetSlicedDataFromROI(data,ROI_low,ROI_high);
        if(i_channel==16 || i_channel==17){continue;}

        int numBins = data.size();

        TH1D* waveformHist = new TH1D("waveformHist", "Waveform", numBins, 0, numBins/5.12);
        for (int i = 0; i < data.size(); ++i)
            waveformHist->Fill(i / 5.12, data[i]/4.095);

        int lowerBin = waveformHist->FindBin(ROI_low/5.12);
        int upperBin = waveformHist->FindBin(ROI_high/5.12);

        float integral = (waveformHist->Integral(lowerBin,upperBin)*50.0*TMath::Power(10.0,-12.0));
        waveformHist->GetXaxis()->SetTitle("Time");
        waveformHist->GetYaxis()->SetTitle("Amplitude");
       

        // Visualize the fit result
        if(vis==1)
        {
            TApplication app("myApp", nullptr, nullptr);
            TCanvas* canvas = new TCanvas("canvas", "Fit Result");
            waveformHist->Draw(); // Draw the original graph
            canvas->Draw();
            app.Run();
        }
        m_data->TD.Charge_Map[i_channel].push_back(integral);

        delete waveformHist;
    }
    return true;
}


bool GetCharge::Finalise()
{
    return true;
}

std::vector<float> GetCharge::GetSlicedDataFromROI(vector<float> data, int startIndex, int endIndex)
{
    if (startIndex < 0 || endIndex > data.size() || startIndex >= endIndex)
    {
        // Invalid range, return an empty vector or handle the error as desired
        return vector<float>();
    }

    // Create a new vector containing the desired slice
    std::vector<float> slice(data.begin() + startIndex, data.begin() + endIndex);

    return slice;
}


void GetCharge::InitRoot()
{
    //m_data->TD.RootFile_Analysis->cd();
    m_data->TD.TTree_Analysis_Charge = new TTree("TTree_Analysis_Charge", "TTree_Analysis_Charge");
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch0", &Charge_Ch0);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch1", &Charge_Ch1);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch2", &Charge_Ch2);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch3", &Charge_Ch3);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch4", &Charge_Ch4);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch5", &Charge_Ch5);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch6", &Charge_Ch6);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch7", &Charge_Ch7);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch8", &Charge_Ch8);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch9", &Charge_Ch9);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch10", &Charge_Ch10);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch11", &Charge_Ch11);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch12", &Charge_Ch12);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch13", &Charge_Ch13);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch14", &Charge_Ch14);
    m_data->TD.TTree_Analysis_Charge->Branch("Charge_Ch15", &Charge_Ch15);        
}