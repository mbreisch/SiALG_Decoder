#include "GetDecayTime.h"

GetDecayTime::GetDecayTime():Tool(){}

bool GetDecayTime::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;


    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
    if(!m_variables.Get("visibility",m_visibility)) m_visibility=0;
    if(!m_variables.Get("ROI_low",ROI_low)) ROI_low=0;
    if(!m_variables.Get("ROI_high",ROI_high)) ROI_high=1023;
    if(!m_variables.Get("MultiplicityCut",MultiplicityCut)) MultiplicityCut=0;

    if(!m_variables.Get("threshold_multiplier",threshold_multiplier)) threshold_multiplier=0.5;      

    if(!m_variables.Get("detector_resolution",m_data->FitModel.detector_resolution)) m_data->FitModel.detector_resolution=0.5;    

    if(!m_variables.Get("guess_tauR",guess_tauR)) guess_tauR=0.5;    
    if(!m_variables.Get("guess_sigma",guess_sigma)) guess_sigma=0.5;    
    if(!m_variables.Get("guess_t0",guess_t0)) guess_t0=0.5;     
    if(!m_variables.Get("guess_n0",guess_n0)) guess_n0=0.5;     
    if(!m_variables.Get("guess_tau0",guess_tau0)) guess_tau0=0.5;     
    if(!m_variables.Get("guess_n1",guess_n1)) guess_n1=0.5;     
    if(!m_variables.Get("guess_tau1",guess_tau1)) guess_tau1=0.5;     
    if(!m_variables.Get("guess_n2",guess_n2)) guess_n2=0.5;     
    if(!m_variables.Get("guess_tau2",guess_tau2)) guess_tau2=0.5;    
    if(!m_variables.Get("guess_n3",guess_n3)) guess_n3=0.5;     
    if(!m_variables.Get("guess_tau3",guess_tau3)) guess_tau3=0.5;  
    if(!m_variables.Get("guess_offset",guess_offset)) guess_offset=0.0;  

    if(!m_variables.Get("guess_gauss_amp",guess_gauss_amp)) guess_gauss_amp=0.5;    
    if(!m_variables.Get("guess_gauss_mean",guess_gauss_mean)) guess_gauss_mean=0.5;     
    if(!m_variables.Get("guess_gauss_sigma",guess_gauss_sigma)) guess_gauss_sigma=0.5;      

    if(!m_variables.Get("fit_type",fit_type)) fit_type="normal";    
    if(!m_variables.Get("fit_options",fit_options)) fit_options="RQ";    
    if(!m_variables.Get("log_style",log_style)) log_style=true;    
    if(!m_variables.Get("normalize",normalize)) normalize=true;    

    if(!m_variables.Get("n_bins",n_bins)) n_bins=1;
    if(!m_variables.Get("bins_start",bins_start)) bins_start=0;    
    if(!m_variables.Get("bins_end",bins_end)) bins_end=1;        

    if(!m_variables.Get("fit_start",fit_start)) fit_start=0;    
    if(!m_variables.Get("fit_end",fit_end)) fit_end=1;        


    return true;
}

bool GetDecayTime::Execute()
{
    if(m_data->TD.Stop==true)
    {
        std::cout <<"Stop has been called, executing Decay calculator"<<std::endl;
        vector<int> temp_channel_list;
        vector<float> CompleteList;

        if(MultiplicityCut==0)
        {
            m_data->TD.Decay_Map = m_data->TD.Decay_Map_uncut;
        }else
        {
            m_data->TD.Decay_Map = m_data->TD.Decay_Map_cut;
        }

        for(const auto& pair : m_data->TD.Decay_Map) 
        {
            temp_channel_list.push_back(pair.first);
        }

        for(int i_channel: temp_channel_list)
        {
            if(i_channel==16 ||i_channel==17){continue;}
            vector<float> c_histogram = m_data->TD.Decay_Map[i_channel];

            CompleteList.insert(CompleteList.end(), c_histogram.begin(), c_histogram.end());

            FitDecay(i_channel,c_histogram);
        }

        FitDecay(404,CompleteList);
        CompleteList.clear();
        return true;
    }

    if(m_data->TD.EndOfRun==true)
    {
        if(MultiplicityCut==0)
        {
            for(const auto& pair : m_data->TD.TTS_Map) 
            {
                int key = pair.first;
                const auto& vector2 = pair.second;

                auto it = m_data->TD.Decay_Map_uncut.find(key);
                if(it != m_data->TD.Decay_Map_uncut.end()) 
                {
                    // Key exists in Decay_Map, extend the vector
                    it->second.insert(it->second.end(), vector2.begin(), vector2.end());
                }else 
                {
                    // Key does not exist in Decay_Map, add a new entry
                    m_data->TD.Decay_Map_uncut[key] = vector2;
                }
            }
        }else
        {
            for(const auto& pair : m_data->TD.TTS_Map_cut) 
            {
                int key = pair.first;
                const auto& vector2 = pair.second;

                auto it = m_data->TD.Decay_Map_cut.find(key);
                if(it != m_data->TD.Decay_Map_cut.end()) 
                {
                    // Key exists in Decay_Map, extend the vector
                    it->second.insert(it->second.end(), vector2.begin(), vector2.end());
                }else 
                {
                    // Key does not exist in Decay_Map, add a new entry
                    m_data->TD.Decay_Map_cut[key] = vector2;
                }
            }
        }
        std::cout <<"Run-End has been called, Decay calculator now has added entries"<<std::endl;
        return true;
    }

    if(m_data->TD.NewRun==true)
    {
        InitRoot();
    }

    return true;
}

bool GetDecayTime::Finalise()
{
    return true;
}

void GetDecayTime::AddLegend(TLegend* legend, const double paramValue, const double paramError, const std::string& paramName) 
{
    std::stringstream stream, streamError;
    stream << std::fixed << std::setprecision(3) << paramValue;
    streamError << std::fixed << std::setprecision(3) << paramError;
    std::string formattedString = paramName + ": " + stream.str() + " +/- " + streamError.str();
    legend->AddEntry((TObject*)0, formattedString.c_str(), "");
} 


///////////////////////////////////////////////////////////////////////
Double_t Sum_of_Exp_Normal_4(Double_t *x, Double_t *par) 
{
    Double_t x_0 = par[0];
    Double_t tau_R = par[1];
    Double_t n_0 = par[2];
    Double_t tau_0 = par[3];
    Double_t n_1 = par[4];
    Double_t tau_1 = par[5];
    Double_t n_2 = par[6];
    Double_t tau_2 = par[7];
    Double_t n_3 = par[8];
    Double_t tau_3 = par[9];

    return n_0 * (  (   TMath::Exp(-(x[0]-x_0)/tau_0) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_0 - tau_R)   ) +
            n_1 * (  (   TMath::Exp(-(x[0]-x_0)/tau_1) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_1 - tau_R)   ) +
            n_2 * (  (   TMath::Exp(-(x[0]-x_0)/tau_2) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_2 - tau_R)   ) +
            n_3 * (  (   TMath::Exp(-(x[0]-x_0)/tau_3) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_3 - tau_R)   ) ;
}
/*
Double_t ExpNew(Double_t *x, Double_t *par) 
{
    Double_t x0 = par[0];
    Double_t n0 = par[1];
    Double_t tau0 = par[2];

    return n0/tau0 * TMath::Exp(-(x[0]-x0)/tau0);
}*/

Double_t ExpNew(Double_t x, Double_t x0, Double_t n0, Double_t tau0) 
{
    return n0 / tau0 * TMath::Exp(-(x - x0) / tau0);
}

Double_t Sum_of_Exp_Fold4(Double_t *x, Double_t *par) 
{
    Double_t x0 = par[0];
    Double_t n_0 = par[1];
    Double_t tau_0 = par[2];
    Double_t n_1 = par[3];
    Double_t tau_1 = par[4];
    Double_t n_2 = par[5];
    Double_t tau_2 = par[6];
    Double_t n_3 = par[7];
    Double_t tau_3 = par[8];

    return (
                (n_0/tau_0) * TMath::Exp(-(x[0] - x0)/tau_0) +
                (n_1/tau_1) * TMath::Exp(-(x[0] - x0)/tau_1) +
                (n_2/tau_2) * TMath::Exp(-(x[0] - x0)/tau_2) + 
                (n_3/tau_3) * TMath::Exp(-(x[0] - x0)/tau_3)
            );
}

Double_t Sum_of_Exp_Normal_3(Double_t *x, Double_t *par) 
{
    Double_t tau_R = par[0];
    Double_t n_0 = par[1];
    Double_t tau_0 = par[2];
    Double_t n_1 = par[3];
    Double_t tau_1 = par[4];
    Double_t n_2 = par[5];
    Double_t tau_2 = par[6];
    Double_t x_0 = par[7];
    Double_t offset = par[8];

    return  n_0 * (  (   TMath::Exp(-(x[0]-x_0)/tau_0) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_0 - tau_R)   ) +
            n_1 * (  (   TMath::Exp(-(x[0]-x_0)/tau_1) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_1 - tau_R)   ) +
            n_2 * (  (   TMath::Exp(-(x[0]-x_0)/tau_2) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_2 - tau_R)   ) + offset;
}

Double_t Sum_of_Exp_Normal_2(Double_t *x, Double_t *par) 
{
    Double_t tau_R = par[0];
    Double_t n_0 = par[1];
    Double_t tau_0 = par[2];
    Double_t n_1 = par[3];
    Double_t tau_1 = par[4];
    Double_t x_0 = par[5];

    return n_0 * (  (   TMath::Exp(-(x[0]-x_0)/tau_0) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_0 - tau_R)   ) +
            n_1 * (  (   TMath::Exp(-(x[0]-x_0)/tau_1) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_1 - tau_R)   );
}

Double_t Single_Exp_Normal(Double_t *x, Double_t *par)
{
    Double_t x_0 = par[0];
    Double_t tau_R = par[1];
    Double_t n_0 = par[2];
    Double_t tau_0 = par[3];
    

    return n_0 * (  (   TMath::Exp(-(x[0]-x_0)/tau_0) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_0 - tau_R)   );
}

Double_t Gauss1(Double_t *x, Double_t *par)
{
    Double_t sigma = par[0];
    return TMath::Exp(-0.5 * TMath::Power(x[0] / sigma, 2));
}

Double_t Gauss_plus_Exp2(Double_t *x, Double_t *par)
{
    Double_t tau_R = par[0];
    Double_t n_0 = par[1];
    Double_t tau_0 = par[2];
    Double_t n_1 = par[3];
    Double_t tau_1 = par[4];
    Double_t x_0 = par[5];
    Double_t offset = par[6];

    Double_t amp = par[7];
    Double_t mean = par[8];
    Double_t sigma = par[9];

    return (amp / (TMath::Sqrt(2.0 * TMath::Pi())*sigma)) * TMath::Exp(-0.5 * TMath::Power((x[0] - mean) / sigma, 2)) +
                (n_0 * (  (   TMath::Exp(-(x[0]-x_0)/tau_0) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_0 - tau_R)   ) +
                 n_1 * (  (   TMath::Exp(-(x[0]-x_0)/tau_1) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_1 - tau_R)   )) + offset;
}

Double_t LinearF0(Double_t *x, Double_t *par)
{
    Double_t offset = par[0];
    return offset;
}

Double_t Sum_of_Exp_Fold(Double_t *x, Double_t *par) 
{
    Double_t sigma = par[0];
    Double_t x0 = par[1];
    Double_t n_0 = par[2];
    Double_t tau_0 = par[3];
    Double_t n_1 = par[4];
    Double_t tau_1 = par[5];
    Double_t n_2 = par[6];
    Double_t tau_2 = par[7];
    Double_t offset = par[8];

    return TMath::Exp(  -TMath::Power(x[0],2) / (2*TMath::Power(sigma,2))   ) * 
            (
                (n_0/tau_0) * TMath::Exp(-(x[0] - x0)/tau_0) +
                (n_1/tau_1) * TMath::Exp(-(x[0] - x0)/tau_1) +
                (n_2/tau_2) * TMath::Exp(-(x[0] - x0)/tau_2)
            ) + offset;
}

Double_t Single_Exp_Fold(Double_t *x, Double_t *par)
{
    Double_t sigma = par[0];
    Double_t x0 = par[1];
    Double_t n_0 = par[2];
    Double_t tau_0 = par[3];

    return  TMath::Exp(  -TMath::Power(x[0],2) / (2*TMath::Power(2.5,2))   ) * 
            (
                (n_0/tau_0) * TMath::Exp(-(x[0] - x0)/tau_0)
            );
}

Double_t Sum_of_Exp_Fold1(Double_t *x, Double_t *par)
{
    Double_t sigma = par[0];
    Double_t x0 = par[1];
    Double_t n_0 = par[2];
    Double_t tau_0 = par[3];
    Double_t offset = par[4];

    return  n_0 * TMath::Exp(-TMath::Power((x[0]-x0),2) / (2*TMath::Power(sigma,2))) *  TMath::Exp(-(x[0])/tau_0);

}

void GetDecayTime::FitDecay(int channel, vector<float> data)
{
    if(data.size()==0){return;}

    TH1F *hist = new TH1F("hist", "Example Histogram", n_bins,bins_start,bins_end);
    for(float k: data) 
    {
        hist->Fill(k);
    }
    if(normalize)
    {
        hist->Scale(1.0 / hist->Integral());
    }
    
    if (fit_type == "gauss") 
    {
        Use_Gaussian(channel, hist);
    }else if (fit_type == "gauss+exp2") 
    {
        Use_Gaussian_plus_Exp2(channel, hist);
    }else if (fit_type == "Exp1") 
    {
        Use_Exp1(channel, hist);
    }else if (fit_type == "Exp2") 
    {
        Use_Exp2(channel, hist);
    }else if (fit_type == "Exp3") 
    {
        Use_Exp3(channel, hist);
    }else if (fit_type == "Exp4") 
    {
        Use_Exp4(channel, hist);
    }else if (fit_type == "Fold1")
    {
        Use_Fold1(channel, hist);
    }else if (fit_type == "Fold3")
    {
        Use_Fold3(channel, hist);
    }else if (fit_type == "Fold4")
    {
        Use_Fold4(channel, hist);
    }else 
    {
        std::cout << RED << "No Fit Type was chosen" << RESET << std::endl;
    }

}

void GetDecayTime::Use_Exp1(int channel, TH1F *hist)
{
    float axis_label_value = (bins_end-bins_start)/n_bins;

    fit_function = new TF1("fit_function", Single_Exp_Normal, fit_start, fit_end, 5); 
    fit_function->SetParameters(guess_t0, guess_n0, guess_tau0, guess_x0);

    // Set bounds for t0 between -100 and 100
    fit_function->SetParLimits(0, 0, 10*guess_t0);
    // Set bounds for n_0 between 0 and 100
    fit_function->SetParLimits(1, 0, 10*guess_n0);
    // Set bounds for tau_0 between 0 and 10
    fit_function->SetParLimits(2, 0, 10*guess_tau0);
    // Set bounds for x_0 between 0 and 10
    fit_function->SetParLimits(3, -100, 100);

    // Perform the fit
    hist->Fit(fit_function, fit_options.c_str()); // "S" for using the sum of squares method

    Double_t fit_t0 = fit_function->GetParameter(0);
    Double_t fit_n_0 = fit_function->GetParameter(1);
    Double_t fit_tau_0 = fit_function->GetParameter(2);
    Double_t fit_x0 = fit_function->GetParameter(3); 

    Double_t fit_t0_error = fit_function->GetParError(0);
    Double_t fit_n_0_error = fit_function->GetParError(1);
    Double_t fit_tau_0_error = fit_function->GetParError(2);
    Double_t fit_x0_error = fit_function->GetParError(3); 


    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << ": t0 = " << fit_t0 << " +/- " << fit_t0_error << " , x0 = " << fit_x0 << " +/- " << fit_x0_error << std::endl;
    std::cout << BLUE << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << RESET << ", " << GREEN << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << RESET << std::endl;
    printf("%s\n","---------------");

    std::string outputName = m_data->TD.Path_Out+"FitData_DecayTime.txt";
    std::ofstream outputFile(outputName.c_str(), std::ios_base::out | std::ios_base::app); // Open a file for writing
    if(outputFile.is_open()) 
    {
        // Write the formatted string to the file
        outputFile << "Fit result, channel " << channel << ": t0 = " << fit_t0 << " +/- " << fit_t0_error << " , x0 = " << fit_x0 << " +/- " << fit_x0_error << std::endl;
        outputFile << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << ", " << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << std::endl;
        outputFile << "---------------" << std::endl;
        outputFile.close(); // Close the file
    }else 
    {
        std::cerr << "Unable to open the file." << std::endl;
    }

    // Create a legend
    TLegend *legend = new TLegend(0.5, 0.6, 0.75, 0.9);
    legend->SetTextSize(0.025);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(2);
    legend->SetMargin(0.1); // Set a 10% margin

    // Add fit parameters and errors to the legend
    std::stringstream t0_stream, t0_stream_error;
    t0_stream << std::fixed << std::setprecision(3) << fit_t0; // Set precision to 3 decimal places
    t0_stream_error << std::fixed << std::setprecision(3) << fit_t0_error; // Set precision to 3 decimal places
    std::string t0_legend = "t0: " + t0_stream.str() + " +/- " + t0_stream_error.str();
    legend->AddEntry((TObject*)0, t0_legend.c_str(), "");

    std::stringstream n0_stream, n0_stream_error;
    n0_stream << std::fixed << std::setprecision(3) << fit_n_0; // Set precision to 3 decimal places
    n0_stream_error << std::fixed << std::setprecision(3) << fit_n_0_error; // Set precision to 3 decimal places
    std::string n0_legend = "n0: " + n0_stream.str() + " +/- " + n0_stream_error.str();
    legend->AddEntry((TObject*)0, n0_legend.c_str(), "");

    std::stringstream tau0_stream, tau0_stream_error;
    tau0_stream << std::fixed << std::setprecision(3) << fit_tau_0; // Set precision to 3 decimal places
    tau0_stream_error << std::fixed << std::setprecision(3) << fit_tau_0_error; // Set precision to 3 decimal places
    std::string tau0_legend = "tau0: " + tau0_stream.str() + " +/- " + tau0_stream_error.str();
    legend->AddEntry((TObject*)0, tau0_legend.c_str(), "");

    std::stringstream label;
    label << std::fixed << std::setprecision(2) << axis_label_value; // Set precision to 3 decimal places

    TCanvas* canvas = new TCanvas("canvas", "Fit Result");
    if(log_style)
    {
        canvas->SetLogy(); 
    }
    hist->SetTitle(("Channel " + std::to_string(channel) + " with multiplicity " + std::to_string(MultiplicityCut)).c_str());
    hist->GetXaxis()->SetTitle("dt between pulse and trigger in ns");
    hist->GetYaxis()->SetTitle((std::string("Counts/")+label.str()+"ns").c_str());
    hist->Draw(); // Draw the original graph
    fit_function->Draw("same"); // Draw the fitted function on top
    legend->Draw();
    canvas->Update();
    std::string name = m_data->TD.Path_Out+"Exp1_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas->SaveAs(name.c_str());  // Save the plot to a file
    canvas->Draw();
    delete canvas;

    // Clean up
    delete hist;
    delete fit_function;    
}

void GetDecayTime::Use_Exp2(int channel, TH1F *hist)
{
    float axis_label_value = (bins_end-bins_start)/n_bins;

    fit_function = new TF1("fit_function", Sum_of_Exp_Normal_2, fit_start, fit_end, 6); 
    fit_function->SetParameters(guess_t0, guess_n0, guess_tau0, guess_n1, guess_tau1, guess_x0);

    // Set bounds for t0 between -100 and 100
    fit_function->SetParLimits(0, 0, 10*guess_t0);
    // Set bounds for n_0 between 0 and 100
    fit_function->SetParLimits(1, 0, 10*guess_n0);
    // Set bounds for tau_0 between 0 and 10
    fit_function->SetParLimits(2, 0, 10*guess_tau0);
    // Set bounds for n_1 between 0 and 100
    fit_function->SetParLimits(3, 0, 10*guess_n1);
    // Set bounds for tau_1 between 0 and 10
    fit_function->SetParLimits(4, 0, 10*guess_tau1);
    // Set bounds for x_0 between 0 and 10
    fit_function->SetParLimits(5, -100, 100);

    // Perform the fit
    hist->Fit(fit_function, fit_options.c_str()); // "S" for using the sum of squares method

    Double_t fit_t0 = fit_function->GetParameter(0);
    Double_t fit_n_0 = fit_function->GetParameter(1);
    Double_t fit_tau_0 = fit_function->GetParameter(2);
    Double_t fit_n_1 = fit_function->GetParameter(3);
    Double_t fit_tau_1 =  fit_function->GetParameter(4);
    Double_t fit_x0 = fit_function->GetParameter(5);

    Double_t fit_t0_error = fit_function->GetParError(0);
    Double_t fit_n_0_error = fit_function->GetParError(1);
    Double_t fit_tau_0_error = fit_function->GetParError(2);
    Double_t fit_n_1_error = fit_function->GetParError(3);
    Double_t fit_tau_1_error = fit_function->GetParError(4);
    Double_t fit_x0_error = fit_function->GetParError(5); 

    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << ": t0 = " << fit_t0 << " +/- " << fit_t0_error << " , x0 = " << fit_x0 << " +/- " << fit_x0_error << std::endl;
    std::cout << BLUE << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << RESET << ", " << GREEN << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << RESET << std::endl;
    std::cout << BLUE << "n1 = " << fit_n_1 << " +/- " << fit_n_1_error << RESET << ", " << GREEN << "tau1 = " << fit_tau_1 << " +/- " << fit_tau_1_error << RESET << std::endl;
    printf("%s\n","---------------");

    std::string outputName = m_data->TD.Path_Out+"FitData_DecayTime.txt";
    std::ofstream outputFile(outputName.c_str(), std::ios_base::out | std::ios_base::app); // Open a file for writing
    if(outputFile.is_open()) 
    {
        // Write the formatted string to the file
        outputFile << "Fit result, channel " << channel << ": t0 = " << fit_t0 << " +/- " << fit_t0_error << " , x0 = " << fit_x0 << " +/- " << fit_x0_error << std::endl;
        outputFile << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << ", " << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << std::endl;
        outputFile << "n1 = " << fit_n_1 << " +/- " << fit_n_1_error << ", " << "tau1 = " << fit_tau_1 << " +/- " << fit_tau_1_error << std::endl;
        outputFile << "---------------" << std::endl;
        outputFile.close(); // Close the file
    }else 
    {
        std::cerr << "Unable to open the file." << std::endl;
    }

    // Create a legend
    TLegend *legend = new TLegend(0.5, 0.6, 0.75, 0.9);
    legend->SetTextSize(0.025);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(2);
    legend->SetMargin(0.1); // Set a 10% margin

    // Add fit parameters and errors to the legend
    std::stringstream t0_stream, t0_stream_error;
    t0_stream << std::fixed << std::setprecision(3) << fit_t0; // Set precision to 3 decimal places
    t0_stream_error << std::fixed << std::setprecision(3) << fit_t0_error; // Set precision to 3 decimal places
    std::string t0_legend = "t0: " + t0_stream.str() + " +/- " + t0_stream_error.str();
    legend->AddEntry((TObject*)0, t0_legend.c_str(), "");

    std::stringstream n0_stream, n0_stream_error;
    n0_stream << std::fixed << std::setprecision(3) << fit_n_0; // Set precision to 3 decimal places
    n0_stream_error << std::fixed << std::setprecision(3) << fit_n_0_error; // Set precision to 3 decimal places
    std::string n0_legend = "n0: " + n0_stream.str() + " +/- " + n0_stream_error.str();
    legend->AddEntry((TObject*)0, n0_legend.c_str(), "");

    std::stringstream tau0_stream, tau0_stream_error;
    tau0_stream << std::fixed << std::setprecision(3) << fit_tau_0; // Set precision to 3 decimal places
    tau0_stream_error << std::fixed << std::setprecision(3) << fit_tau_0_error; // Set precision to 3 decimal places
    std::string tau0_legend = "tau0: " + tau0_stream.str() + " +/- " + tau0_stream_error.str();
    legend->AddEntry((TObject*)0, tau0_legend.c_str(), "");

    std::stringstream n1_stream, n1_stream_error;
    n1_stream << std::fixed << std::setprecision(3) << fit_n_1; // Set precision to 3 decimal places
    n1_stream_error << std::fixed << std::setprecision(3) << fit_n_1_error; // Set precision to 3 decimal places
    std::string n1_legend = "n1: " + n1_stream.str() + " +/- " + n1_stream_error.str();
    legend->AddEntry((TObject*)0, n1_legend.c_str(), "");

    std::stringstream tau1_stream, tau1_stream_error;
    tau1_stream << std::fixed << std::setprecision(3) << fit_tau_1; // Set precision to 3 decimal places
    tau1_stream_error << std::fixed << std::setprecision(3) << fit_tau_1_error; // Set precision to 3 decimal places
    std::string tau1_legend = "tau1: " + tau1_stream.str() + " +/- " + tau1_stream_error.str();
    legend->AddEntry((TObject*)0, tau1_legend.c_str(), "");

    std::stringstream label;
    label << std::fixed << std::setprecision(2) << axis_label_value; // Set precision to 3 decimal places

    singleExpFunc0 = new TF1("singleExpFunc0", Single_Exp_Normal, fit_start, fit_end, 4); // 4 parameters for Single_Exp
    singleExpFunc1 = new TF1("singleExpFunc1", Single_Exp_Normal, fit_start, fit_end, 4); // 4 parameters for Single_Exp

    singleExpFunc0->SetParameters(fit_t0, fit_n_0, fit_tau_0, fit_x0);
    singleExpFunc1->SetParameters(fit_t0, fit_n_1, fit_tau_1, fit_x0);

    singleExpFunc0->SetLineColor(kGray);
    singleExpFunc0->SetLineStyle(2);
    singleExpFunc1->SetLineColor(kGreen);
    singleExpFunc1->SetLineStyle(2);

    TCanvas* canvas = new TCanvas("canvas", "Fit Result");
    if(log_style)
    {
        canvas->SetLogy(); 
    }
    hist->SetTitle(("Channel " + std::to_string(channel) + " with multiplicity " + std::to_string(MultiplicityCut)).c_str());
    hist->GetXaxis()->SetTitle("dt between pulse and trigger in ns");
    hist->GetYaxis()->SetTitle((std::string("Counts/")+label.str()+"ns").c_str());
    hist->Draw(); // Draw the original graph
    fit_function->Draw("same"); // Draw the fitted function on top
    singleExpFunc0->Draw("same");
    singleExpFunc1->Draw("same");
    legend->Draw();
    canvas->Update();
    std::string name = m_data->TD.Path_Out+"Exp2_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas->SaveAs(name.c_str());  // Save the plot to a file
    canvas->Draw();
    delete canvas;

    // Clean up
    delete hist;
    delete fit_function;
}

void GetDecayTime::Use_Exp3(int channel, TH1F *hist)
{
    float axis_label_value = (bins_end-bins_start)/n_bins;

    fit_function = new TF1("fit_function", Sum_of_Exp_Normal_3, fit_start, fit_end, 9); 
    fit_function->SetParameters(guess_t0, guess_n0, guess_tau0, guess_n1, guess_tau1, guess_n2, guess_tau2, guess_x0, guess_offset);

    // Set bounds for t0 between -100 and 100
    fit_function->SetParLimits(0, 0, 10*guess_t0);
    // Set bounds for n_0 between 0 and 100
    fit_function->SetParLimits(1, 0, 10*guess_n0);
    // Set bounds for tau_0 between 0 and 10
    fit_function->SetParLimits(2, 0, 10*guess_tau0);
    // Set bounds for n_1 between 0 and 100
    fit_function->SetParLimits(3, 0, 10*guess_n1);
    // Set bounds for tau_1 between 0 and 10
    fit_function->SetParLimits(4, 0, 10*guess_tau1);
    // Set bounds for n_2 between 0 and 100
    fit_function->SetParLimits(5, 0, 10*guess_n2);
    // Set bounds for tau_2 between 0 and 10
    fit_function->SetParLimits(6, 0, 10*guess_tau2);
    // Set bounds for x_0 between 0 and 10
    fit_function->SetParLimits(7, -100, 100);
    // Set bounds for offset between 0 and 10
    fit_function->SetParLimits(8, 0, 10*guess_offset);

    // Perform the fit
    hist->Fit(fit_function, fit_options.c_str()); // "S" for using the sum of squares method

    Double_t fit_t0 = fit_function->GetParameter(0);
    Double_t fit_n_0 = fit_function->GetParameter(1);
    Double_t fit_tau_0 = fit_function->GetParameter(2);
    Double_t fit_n_1 = fit_function->GetParameter(3);
    Double_t fit_tau_1 =  fit_function->GetParameter(4);
    Double_t fit_n_2 = fit_function->GetParameter(5);
    Double_t fit_tau_2 = fit_function->GetParameter(6);
    Double_t fit_x0 = fit_function->GetParameter(7); 
    Double_t fit_offset = fit_function->GetParameter(8); 

    Double_t fit_t0_error = fit_function->GetParError(0);
    Double_t fit_n_0_error = fit_function->GetParError(1);
    Double_t fit_tau_0_error = fit_function->GetParError(2);
    Double_t fit_n_1_error = fit_function->GetParError(3);
    Double_t fit_tau_1_error = fit_function->GetParError(4);
    Double_t fit_n_2_error = fit_function->GetParError(5);
    Double_t fit_tau_2_error = fit_function->GetParError(6);
    Double_t fit_x0_error = fit_function->GetParError(7); 

    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << ": t0 = " << fit_t0 << " +/- " << fit_t0_error << " , x0 = " << fit_x0 << " +/- " << fit_x0_error << std::endl;
    std::cout << BLUE << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << RESET << ", " << GREEN << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << RESET << std::endl;
    std::cout << BLUE << "n1 = " << fit_n_1 << " +/- " << fit_n_1_error << RESET << ", " << GREEN << "tau1 = " << fit_tau_1 << " +/- " << fit_tau_1_error << RESET << std::endl;
    std::cout << BLUE << "n2 = " << fit_n_2 << " +/- " << fit_n_2_error << RESET << ", " << GREEN << "tau2 = " << fit_tau_2 << " +/- " << fit_tau_2_error << RESET << std::endl;
    printf("%s\n","---------------");

    std::string outputName = m_data->TD.Path_Out+"FitData_DecayTime.txt";
    std::ofstream outputFile(outputName.c_str(), std::ios_base::out | std::ios_base::app); // Open a file for writing
    if(outputFile.is_open()) 
    {
        // Write the formatted string to the file
        outputFile << "Fit result, channel " << channel << ": t0 = " << fit_t0 << " +/- " << fit_t0_error << " , x0 = " << fit_x0 << " +/- " << fit_x0_error << std::endl;
        outputFile << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << ", " << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << std::endl;
        outputFile << "n1 = " << fit_n_1 << " +/- " << fit_n_1_error << ", " << "tau1 = " << fit_tau_1 << " +/- " << fit_tau_1_error << std::endl;
        outputFile << "n2 = " << fit_n_2 << " +/- " << fit_n_2_error << ", " << "tau2 = " << fit_tau_2 << " +/- " << fit_tau_2_error << std::endl;
        outputFile << "---------------" << std::endl;
        outputFile.close(); // Close the file
    }else 
    {
        std::cerr << "Unable to open the file." << std::endl;
    }

    // Create a legend
    TLegend *legend = new TLegend(0.5, 0.6, 0.75, 0.9);
    legend->SetTextSize(0.025);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(2);
    legend->SetMargin(0.1); // Set a 10% margin

    // Add fit parameters and errors to the legend
    std::stringstream t0_stream, t0_stream_error;
    t0_stream << std::fixed << std::setprecision(3) << fit_t0; // Set precision to 3 decimal places
    t0_stream_error << std::fixed << std::setprecision(3) << fit_t0_error; // Set precision to 3 decimal places
    std::string t0_legend = "t0: " + t0_stream.str() + " +/- " + t0_stream_error.str();
    legend->AddEntry((TObject*)0, t0_legend.c_str(), "");

    std::stringstream n0_stream, n0_stream_error;
    n0_stream << std::fixed << std::setprecision(3) << fit_n_0; // Set precision to 3 decimal places
    n0_stream_error << std::fixed << std::setprecision(3) << fit_n_0_error; // Set precision to 3 decimal places
    std::string n0_legend = "n0: " + n0_stream.str() + " +/- " + n0_stream_error.str();
    legend->AddEntry((TObject*)0, n0_legend.c_str(), "");

    std::stringstream tau0_stream, tau0_stream_error;
    tau0_stream << std::fixed << std::setprecision(3) << fit_tau_0; // Set precision to 3 decimal places
    tau0_stream_error << std::fixed << std::setprecision(3) << fit_tau_0_error; // Set precision to 3 decimal places
    std::string tau0_legend = "tau0: " + tau0_stream.str() + " +/- " + tau0_stream_error.str();
    legend->AddEntry((TObject*)0, tau0_legend.c_str(), "");

    std::stringstream n1_stream, n1_stream_error;
    n1_stream << std::fixed << std::setprecision(3) << fit_n_1; // Set precision to 3 decimal places
    n1_stream_error << std::fixed << std::setprecision(3) << fit_n_1_error; // Set precision to 3 decimal places
    std::string n1_legend = "n1: " + n1_stream.str() + " +/- " + n1_stream_error.str();
    legend->AddEntry((TObject*)0, n1_legend.c_str(), "");

    std::stringstream tau1_stream, tau1_stream_error;
    tau1_stream << std::fixed << std::setprecision(3) << fit_tau_1; // Set precision to 3 decimal places
    tau1_stream_error << std::fixed << std::setprecision(3) << fit_tau_1_error; // Set precision to 3 decimal places
    std::string tau1_legend = "tau1: " + tau1_stream.str() + " +/- " + tau1_stream_error.str();
    legend->AddEntry((TObject*)0, tau1_legend.c_str(), "");

    std::stringstream n2_stream, n2_stream_error;
    n2_stream << std::fixed << std::setprecision(3) << fit_n_2; // Set precision to 3 decimal places
    n2_stream_error << std::fixed << std::setprecision(3) << fit_n_2_error; // Set precision to 3 decimal places
    std::string n2_legend = "n2: " + n2_stream.str() + " +/- " + n2_stream_error.str();
    legend->AddEntry((TObject*)0, n2_legend.c_str(), "");

    std::stringstream tau2_stream, tau2_stream_error;
    tau2_stream << std::fixed << std::setprecision(3) << fit_tau_2; // Set precision to 3 decimal places
    tau2_stream_error << std::fixed << std::setprecision(3) << fit_tau_2_error; // Set precision to 3 decimal places
    std::string tau2_legend = "tau2: " + tau2_stream.str() + " +/- " + tau2_stream_error.str();
    legend->AddEntry((TObject*)0, tau2_legend.c_str(), "");

    std::stringstream label;
    label << std::fixed << std::setprecision(2) << axis_label_value; // Set precision to 3 decimal places

    singleExpFunc0 = new TF1("singleExpFunc0", Single_Exp_Normal, bins_start, bins_end, 4); // 4 parameters for Single_Exp
    singleExpFunc1 = new TF1("singleExpFunc1", Single_Exp_Normal, bins_start, bins_end, 4); // 4 parameters for Single_Exp
    singleExpFunc2 = new TF1("singleExpFunc2", Single_Exp_Normal, bins_start, bins_end, 4); // 4 parameters for Single_Exp
    linearFunc0 = new TF1("linearFunc0", LinearF0, bins_start, bins_end, 5); // 4 parameters for Single_Exp

    singleExpFunc0->SetParameters(fit_t0, fit_n_0, fit_tau_0, fit_x0);
    singleExpFunc1->SetParameters(fit_t0, fit_n_1, fit_tau_1, fit_x0);
    singleExpFunc2->SetParameters(fit_t0, fit_n_2, fit_tau_2, fit_x0);
    linearFunc0->SetParameter(0,fit_offset);

    singleExpFunc0->SetLineColor(kGray);
    singleExpFunc0->SetLineStyle(2);
    singleExpFunc1->SetLineColor(kGreen);
    singleExpFunc1->SetLineStyle(2);
    singleExpFunc2->SetLineColor(kViolet);
    singleExpFunc2->SetLineStyle(2);
    linearFunc0->SetLineColor(kYellow);
    linearFunc0->SetLineStyle(2);

    TCanvas* canvas = new TCanvas("canvas", "Fit Result");
    if(log_style)
    {
        canvas->SetLogy(); 
    }
    hist->SetTitle(("Channel " + std::to_string(channel) + " with multiplicity " + std::to_string(MultiplicityCut)).c_str());
    hist->GetXaxis()->SetTitle("dt between pulse and trigger in ns");
    hist->GetYaxis()->SetTitle((std::string("Counts/")+label.str()+"ns").c_str());
    hist->Draw(); // Draw the original graph
    fit_function->Draw("same"); // Draw the fitted function on top
    singleExpFunc0->Draw("same");
    singleExpFunc1->Draw("same");
    singleExpFunc2->Draw("same");
    linearFunc0->Draw("same");
    legend->Draw();
    canvas->Update();
    std::string name = m_data->TD.Path_Out+"Exp3_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas->SaveAs(name.c_str());  // Save the plot to a file
    canvas->Draw();
    delete canvas;

    // Clean up
    delete hist;
    delete fit_function;
}

void GetDecayTime::Use_Exp4(int channel, TH1F *hist)
{
    float axis_label_value = (bins_end-bins_start)/n_bins;
    m_data->FitModel.summation_order = 4;

    // Create a workspace
    RooWorkspace w("w");

    // Define the observable variable
    RooRealVar x("x", "x", -10, 10);

    // Create custom PDFs using RooGenericPdf
    RooRealVar mean("mean", "mean", 0);  // Set your desired mean value here
    RooRealVar sigma("sigma", "sigma", 1.5); 
    RooGaussian gauss("gauss", "gauss", x, mean, sigma);

    RooRealVar par0("par0", "par0", guess_t0);
    RooRealVar par1("par1", "par1", guess_n0);
    RooRealVar par2("par2", "par2", guess_tau0);
    RooArgList pdfComponents;
    pdfComponents.add(RooArgSet(x));
    pdfComponents.add(par0);
    pdfComponents.add(par1);
    pdfComponents.add(par2);

    // Creating the RooGenericPdf
    RooGenericPdf expNew("expNew", "ExpNew", "par1/par2 * exp(-(x-par0)/par2)", pdfComponents);

    RooFFTConvPdf convolution("convolution", "convolution", x, gauss, expNew);

    // Convert the TH1F histogram to a RooDataHist
    RooDataHist dataHist("dataHist", "Data Hist", RooArgSet(x), hist);

    // Perform the fit to the RooDataHist
    //convolution.fitTo(dataHist);

    // Plot the result
    RooPlot* frame = x.frame();
    dataHist.plotOn(frame);
    gauss.plotOn(frame, RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
    //expNew.plotOn(frame, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
    //convolution.plotOn(frame, RooFit::LineColor(kGreen));

    // Create a canvas and draw the frame
    TCanvas canvas("canvas", "Convolution Fit Example", 800, 600);
    frame->Draw();
    std::string name = m_data->TD.Path_Out+"Conv_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas.SaveAs(name.c_str());

    std::cout<<"---------------------------"<<std::endl;
}

void GetDecayTime::Use_Fold3(int channel, TH1F *hist)
{
    float axis_label_value = (bins_end-bins_start)/n_bins;

    fit_function = new TF1("fit_function", Sum_of_Exp_Fold, fit_start, fit_end, 9); 
    fit_function->SetParameters(guess_sigma, guess_x0, guess_n0, guess_tau0, guess_n1, guess_tau1, guess_n2, guess_tau2, guess_offset);
    
    // Set bounds for sigma between -100 and 100
    fit_function->SetParLimits(0, 0, 10*guess_sigma);
    // Set bounds for x0 between -100 and 100
    fit_function->SetParLimits(1, 10*guess_x0, 0);
    // Set bounds for n_0 between 0 and 100
    fit_function->SetParLimits(2, 0, 10*guess_n0);
    // Set bounds for tau_0 between 0 and 10
    fit_function->SetParLimits(3, 0, 10*guess_tau0);
    // Set bounds for n_1 between 0 and 100
    fit_function->SetParLimits(4, 0, 10*guess_n1);
    // Set bounds for tau_1 between 0 and 10
    fit_function->SetParLimits(5, 0, 10*guess_tau1);
    // Set bounds for n_2 between 0 and 100
    fit_function->SetParLimits(6, 0, 10*guess_n2);
    // Set bounds for tau_2 between 0 and 10
    fit_function->SetParLimits(7, 0, 10*guess_tau2);
    // Set bounds for offset between 0 and 10
    fit_function->SetParLimits(8, 0, 10*guess_offset);

    // Perform the fit
    hist->Fit(fit_function, fit_options.c_str()); // "S" for using the sum of squares method

    Double_t fit_sigma = fit_function->GetParameter(0);
    Double_t fit_x0 = fit_function->GetParameter(1); 
    Double_t fit_n_0 = fit_function->GetParameter(2);
    Double_t fit_tau_0 = fit_function->GetParameter(3);
    Double_t fit_n_1 = fit_function->GetParameter(4);
    Double_t fit_tau_1 =  fit_function->GetParameter(5);
    Double_t fit_n_2 = fit_function->GetParameter(6);
    Double_t fit_tau_2 = fit_function->GetParameter(7);
    Double_t fit_offset = fit_function->GetParameter(8); 

    Double_t fit_sigma_error = fit_function->GetParError(0);
    Double_t fit_x0_error = fit_function->GetParError(1); 
    Double_t fit_n_0_error = fit_function->GetParError(2);
    Double_t fit_tau_0_error = fit_function->GetParError(3);
    Double_t fit_n_1_error = fit_function->GetParError(4);
    Double_t fit_tau_1_error = fit_function->GetParError(5);
    Double_t fit_n_2_error = fit_function->GetParError(6);
    Double_t fit_tau_2_error = fit_function->GetParError(7);
    Double_t fit_offset_error = fit_function->GetParameter(8);
    

    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << ": x0 = " << fit_x0 << " +/- " << fit_x0_error << " , off = " << fit_offset << " +/- " << fit_offset_error << std::endl;
    std::cout << BLUE << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << RESET << ", " << GREEN << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << RESET << std::endl;
    std::cout << BLUE << "n1 = " << fit_n_1 << " +/- " << fit_n_1_error << RESET << ", " << GREEN << "tau1 = " << fit_tau_1 << " +/- " << fit_tau_1_error << RESET << std::endl;
    std::cout << BLUE << "n2 = " << fit_n_2 << " +/- " << fit_n_2_error << RESET << ", " << GREEN << "tau2 = " << fit_tau_2 << " +/- " << fit_tau_2_error << RESET << std::endl;
    printf("%s\n","---------------");

    std::string outputName = m_data->TD.Path_Out+"FitData_DecayTime.txt";
    std::ofstream outputFile(outputName.c_str(), std::ios_base::out | std::ios_base::app); // Open a file for writing
    if(outputFile.is_open()) 
    {
        // Write the formatted string to the file
        outputFile << "Fit result, channel " << channel << ": x0 = " << fit_x0 << " +/- " << fit_x0_error << " , off = " << fit_offset << " +/- " << fit_offset_error << std::endl;
        outputFile << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << ", " << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << std::endl;
        outputFile << "n1 = " << fit_n_1 << " +/- " << fit_n_1_error << ", " << "tau1 = " << fit_tau_1 << " +/- " << fit_tau_1_error << std::endl;
        outputFile << "n2 = " << fit_n_2 << " +/- " << fit_n_2_error << ", " << "tau2 = " << fit_tau_2 << " +/- " << fit_tau_2_error << std::endl;
        outputFile << "---------------" << std::endl;
        outputFile.close(); // Close the file
    }else 
    {
        std::cerr << "Unable to open the file." << std::endl;
    }

    // Create a legend
    TLegend *legend = new TLegend(0.5, 0.6, 0.75, 0.9);
    legend->SetTextSize(0.025);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(2);
    legend->SetMargin(0.1); // Set a 10% margin

    // Add fit parameters and errors to the legend
    std::stringstream sigma_stream, sigma_stream_error;
    sigma_stream << std::fixed << std::setprecision(3) << fit_sigma; // Set precision to 3 decimal places
    sigma_stream_error << std::fixed << std::setprecision(3) << fit_sigma_error; // Set precision to 3 decimal places
    std::string sigma_legend = "sigma: " + sigma_stream.str() + " +/- " + sigma_stream_error.str();
    legend->AddEntry((TObject*)0, sigma_legend.c_str(), "");

    std::stringstream x0_stream, x0_stream_error;
    x0_stream << std::fixed << std::setprecision(3) << fit_x0; // Set precision to 3 decimal places
    x0_stream_error << std::fixed << std::setprecision(3) << fit_x0_error; // Set precision to 3 decimal places
    std::string x0_legend = "x0: " + x0_stream.str() + " +/- " + x0_stream_error.str();
    legend->AddEntry((TObject*)0, x0_legend.c_str(), "");

    std::stringstream n0_stream, n0_stream_error;
    n0_stream << std::fixed << std::setprecision(3) << fit_n_0; // Set precision to 3 decimal places
    n0_stream_error << std::fixed << std::setprecision(3) << fit_n_0_error; // Set precision to 3 decimal places
    std::string n0_legend = "n0: " + n0_stream.str() + " +/- " + n0_stream_error.str();
    legend->AddEntry((TObject*)0, n0_legend.c_str(), "");

    std::stringstream tau0_stream, tau0_stream_error;
    tau0_stream << std::fixed << std::setprecision(3) << fit_tau_0; // Set precision to 3 decimal places
    tau0_stream_error << std::fixed << std::setprecision(3) << fit_tau_0_error; // Set precision to 3 decimal places
    std::string tau0_legend = "tau0: " + tau0_stream.str() + " +/- " + tau0_stream_error.str();
    legend->AddEntry((TObject*)0, tau0_legend.c_str(), "");

    std::stringstream n1_stream, n1_stream_error;
    n1_stream << std::fixed << std::setprecision(3) << fit_n_1; // Set precision to 3 decimal places
    n1_stream_error << std::fixed << std::setprecision(3) << fit_n_1_error; // Set precision to 3 decimal places
    std::string n1_legend = "n1: " + n1_stream.str() + " +/- " + n1_stream_error.str();
    legend->AddEntry((TObject*)0, n1_legend.c_str(), "");

    std::stringstream tau1_stream, tau1_stream_error;
    tau1_stream << std::fixed << std::setprecision(3) << fit_tau_1; // Set precision to 3 decimal places
    tau1_stream_error << std::fixed << std::setprecision(3) << fit_tau_1_error; // Set precision to 3 decimal places
    std::string tau1_legend = "tau1: " + tau1_stream.str() + " +/- " + tau1_stream_error.str();
    legend->AddEntry((TObject*)0, tau1_legend.c_str(), "");

    std::stringstream n2_stream, n2_stream_error;
    n2_stream << std::fixed << std::setprecision(3) << fit_n_2; // Set precision to 3 decimal places
    n2_stream_error << std::fixed << std::setprecision(3) << fit_n_2_error; // Set precision to 3 decimal places
    std::string n2_legend = "n2: " + n2_stream.str() + " +/- " + n2_stream_error.str();
    legend->AddEntry((TObject*)0, n2_legend.c_str(), "");

    std::stringstream tau2_stream, tau2_stream_error;
    tau2_stream << std::fixed << std::setprecision(3) << fit_tau_2; // Set precision to 3 decimal places
    tau2_stream_error << std::fixed << std::setprecision(3) << fit_tau_2_error; // Set precision to 3 decimal places
    std::string tau2_legend = "tau2: " + tau2_stream.str() + " +/- " + tau2_stream_error.str();
    legend->AddEntry((TObject*)0, tau2_legend.c_str(), "");

    std::stringstream offset_stream, offset_stream_error;
    offset_stream << std::fixed << std::setprecision(3) << fit_offset; // Set precision to 3 decimal places
    offset_stream_error << std::fixed << std::setprecision(3) << fit_offset_error; // Set precision to 3 decimal places
    std::string offset_legend = "offset: " + offset_stream.str() + " +/- " + offset_stream_error.str();
    legend->AddEntry((TObject*)0, offset_legend.c_str(), "");

    std::stringstream label;
    label << std::fixed << std::setprecision(2) << axis_label_value; // Set precision to 3 decimal places

    singleExpFunc0 = new TF1("singleExpFunc0", Single_Exp_Fold, bins_start, bins_end, 4); // 4 parameters for Single_Exp
    singleExpFunc1 = new TF1("singleExpFunc1", Single_Exp_Fold, bins_start, bins_end, 4); // 4 parameters for Single_Exp
    singleExpFunc2 = new TF1("singleExpFunc2", Single_Exp_Fold, bins_start, bins_end, 4); // 4 parameters for Single_Exp
    linearFunc0 = new TF1("linearFunc0", LinearF0, bins_start, bins_end, 1); // 4 parameters for Single_Exp

    singleExpFunc0->SetParameters(fit_sigma, fit_x0, fit_n_0, fit_tau_0);
    singleExpFunc1->SetParameters(fit_sigma, fit_x0, fit_n_1, fit_tau_1);
    singleExpFunc2->SetParameters(fit_sigma, fit_x0, fit_n_2, fit_tau_2);
    linearFunc0->SetParameter(0,fit_offset);

    singleExpFunc0->SetLineColor(kGray);
    singleExpFunc0->SetLineStyle(2);
    singleExpFunc1->SetLineColor(kGreen);
    singleExpFunc1->SetLineStyle(2);
    singleExpFunc2->SetLineColor(kViolet);
    singleExpFunc2->SetLineStyle(2);
    linearFunc0->SetLineColor(kYellow);
    linearFunc0->SetLineStyle(2);

    TCanvas* canvas = new TCanvas("canvas", "Fit Result", 800, 600);
    if(log_style)
    {
        canvas->SetLogy(); 
    }
    hist->SetTitle(("Channel " + std::to_string(channel) + " with multiplicity " + std::to_string(MultiplicityCut)).c_str());
    hist->GetXaxis()->SetTitle("dt between pulse and trigger in ns");
    hist->GetYaxis()->SetTitle((std::string("Counts/")+label.str()+"ns").c_str());
    hist->Draw(); // Draw the original graph
    fit_function->Draw("same"); // Draw the fitted function on top
    singleExpFunc0->Draw("same");
    singleExpFunc1->Draw("same");
    singleExpFunc2->Draw("same");
    linearFunc0->Draw("same");
    legend->Draw();
    canvas->Update();
    std::string name = m_data->TD.Path_Out+"Fold3_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas->SaveAs(name.c_str());  // Save the plot to a file
    canvas->Draw();
    delete canvas;

    // Clean up
    delete hist;
    delete fit_function;
}

void GetDecayTime::Use_Fold1(int channel, TH1F *hist)
{
    float axis_label_value = (bins_end-bins_start)/n_bins;

    fit_function = new TF1("fit_function", Sum_of_Exp_Fold1, fit_start, fit_end, 5); 
    fit_function->SetParameters(guess_sigma, guess_x0, guess_n0, guess_tau0, guess_offset);
    
    // Set bounds for sigma between -100 and 100
    fit_function->SetParLimits(0, 0, 10*guess_sigma);
    // Set bounds for x0 between -100 and 100
    fit_function->SetParLimits(1, 0, 10*guess_x0);
    // Set bounds for n_0 between 0 and 100
    fit_function->SetParLimits(2, 0, 10*guess_n0);
    // Set bounds for tau_0 between 0 and 10
    fit_function->SetParLimits(3, 0, 10*guess_tau0);
    // Set bounds for offset between 0 and 10
    fit_function->SetParLimits(4, 0, 10*guess_offset);

    // Perform the fit
    hist->Fit(fit_function, fit_options.c_str()); // "S" for using the sum of squares method

    Double_t fit_sigma = fit_function->GetParameter(0);
    Double_t fit_x0 = fit_function->GetParameter(1); 
    Double_t fit_n_0 = fit_function->GetParameter(2);
    Double_t fit_tau_0 = fit_function->GetParameter(3);
    Double_t fit_offset = fit_function->GetParameter(4); 

    Double_t fit_sigma_error = fit_function->GetParError(0);
    Double_t fit_x0_error = fit_function->GetParError(1); 
    Double_t fit_n_0_error = fit_function->GetParError(2);
    Double_t fit_tau_0_error = fit_function->GetParError(3);
    Double_t fit_offset_error = fit_function->GetParameter(4);
    
    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << ": x0 = " << fit_x0 << " +/- " << fit_x0_error << " , off = " << fit_offset << " +/- " << fit_offset_error << std::endl;
    std::cout << BLUE << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << RESET << ", " << GREEN << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << RESET << std::endl;
    printf("%s\n","---------------");

    std::string outputName = m_data->TD.Path_Out+"FitData_DecayTime.txt";
    std::ofstream outputFile(outputName.c_str(), std::ios_base::out | std::ios_base::app); // Open a file for writing
    if(outputFile.is_open()) 
    {
        // Write the formatted string to the file
        outputFile << "Fit result, channel " << channel << ": x0 = " << fit_x0 << " +/- " << fit_x0_error << " , off = " << fit_offset << " +/- " << fit_offset_error << std::endl;
        outputFile << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << ", " << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << std::endl;
        outputFile << "---------------" << std::endl;
        outputFile.close(); // Close the file
    }else 
    {
        std::cerr << "Unable to open the file." << std::endl;
    }

    // Create a legend
    TLegend *legend = new TLegend(0.5, 0.6, 0.75, 0.9);
    legend->SetTextSize(0.025);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(2);
    legend->SetMargin(0.1); // Set a 10% margin

    // Add fit parameters and errors to the legend
    std::stringstream sigma_stream, sigma_stream_error;
    sigma_stream << std::fixed << std::setprecision(3) << fit_sigma; // Set precision to 3 decimal places
    sigma_stream_error << std::fixed << std::setprecision(3) << fit_sigma_error; // Set precision to 3 decimal places
    std::string sigma_legend = "sigma: " + sigma_stream.str() + " +/- " + sigma_stream_error.str();
    legend->AddEntry((TObject*)0, sigma_legend.c_str(), "");

    std::stringstream x0_stream, x0_stream_error;
    x0_stream << std::fixed << std::setprecision(3) << fit_x0; // Set precision to 3 decimal places
    x0_stream_error << std::fixed << std::setprecision(3) << fit_x0_error; // Set precision to 3 decimal places
    std::string x0_legend = "x0: " + x0_stream.str() + " +/- " + x0_stream_error.str();
    legend->AddEntry((TObject*)0, x0_legend.c_str(), "");

    std::stringstream n0_stream, n0_stream_error;
    n0_stream << std::fixed << std::setprecision(3) << fit_n_0; // Set precision to 3 decimal places
    n0_stream_error << std::fixed << std::setprecision(3) << fit_n_0_error; // Set precision to 3 decimal places
    std::string n0_legend = "n0: " + n0_stream.str() + " +/- " + n0_stream_error.str();
    legend->AddEntry((TObject*)0, n0_legend.c_str(), "");

    std::stringstream tau0_stream, tau0_stream_error;
    tau0_stream << std::fixed << std::setprecision(3) << fit_tau_0; // Set precision to 3 decimal places
    tau0_stream_error << std::fixed << std::setprecision(3) << fit_tau_0_error; // Set precision to 3 decimal places
    std::string tau0_legend = "tau0: " + tau0_stream.str() + " +/- " + tau0_stream_error.str();
    legend->AddEntry((TObject*)0, tau0_legend.c_str(), "");

    std::stringstream offset_stream, offset_stream_error;
    offset_stream << std::fixed << std::setprecision(3) << fit_offset; // Set precision to 3 decimal places
    offset_stream_error << std::fixed << std::setprecision(3) << fit_offset_error; // Set precision to 3 decimal places
    std::string offset_legend = "offset: " + offset_stream.str() + " +/- " + offset_stream_error.str();
    legend->AddEntry((TObject*)0, offset_legend.c_str(), "");

    std::stringstream label;
    label << std::fixed << std::setprecision(2) << axis_label_value; // Set precision to 3 decimal places

    linearFunc0 = new TF1("linearFunc0", LinearF0, bins_start, bins_end, 1); // 4 parameters for Single_Exp
    linearFunc0->SetParameter(0,fit_offset);
    linearFunc0->SetLineColor(kYellow);
    linearFunc0->SetLineStyle(2);

    TCanvas* canvas = new TCanvas("canvas", "Fit Result", 800, 600);
    if(log_style)
    {
        canvas->SetLogy(); 
    }
    hist->SetTitle(("Channel " + std::to_string(channel) + " with multiplicity " + std::to_string(MultiplicityCut)).c_str());
    hist->GetXaxis()->SetTitle("dt between pulse and trigger in ns");
    hist->GetYaxis()->SetTitle((std::string("Counts/")+label.str()+"ns").c_str());
    hist->Draw(); // Draw the original graph
    fit_function->Draw("same"); // Draw the fitted function on top
    linearFunc0->Draw("same");
    legend->Draw();
    canvas->Update();
    std::string name = m_data->TD.Path_Out+"Fold1_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas->SaveAs(name.c_str());  // Save the plot to a file
    canvas->Draw();
    delete canvas;

    // Clean up
    delete hist;
    delete fit_function;
}

void GetDecayTime::Use_Gaussian(int channel, TH1F *hist)
{
    int axis_label_value = (bins_end-bins_start)/n_bins;

    fit_function = new TF1("Gaussian", FitModels::Gaussian, fit_start, fit_end, 3);
    fit_function->SetParameters(guess_gauss_amp,guess_gauss_mean,guess_gauss_sigma);

    fit_function->SetParLimits(0,   0,      1);
    fit_function->SetParLimits(1,   -100,   100);
    fit_function->SetParLimits(2,   0,      10*guess_gauss_sigma);

    hist->Fit(fit_function, fit_options.c_str());

    Double_t Fit_Amp = fit_function->GetParameter(0);
    Double_t Fit_Amp_e = fit_function->GetParError(0);

    Double_t Fit_Mean = fit_function->GetParameter(1);
    Double_t Fit_Mean_e = fit_function->GetParError(1);

    Double_t Fit_Sigma = fit_function->GetParameter(2);
    Double_t Fit_Sigma_e = fit_function->GetParError(2);

    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << std::endl;
    std::cout << BLUE << "Amp = " << Fit_Amp << " +/- " << Fit_Amp_e << RESET << std::endl;
    std::cout << BLUE << "Mean = " << Fit_Mean << " +/- " << Fit_Mean_e << RESET << std::endl;
    std::cout << GREEN << "Sigma = " << Fit_Sigma << " +/- " << Fit_Sigma_e << RESET << std::endl;
    printf("%s\n","---------------");

    std::string outputName = m_data->TD.Path_Out+"FitData_DecayTime.txt";
    std::ofstream outputFile(outputName.c_str(), std::ios_base::out | std::ios_base::app); // Open a file for writing
    if(outputFile.is_open()) 
    {
        // Write the formatted string to the file
        outputFile << "Fit result, channel " << channel << ", fit type was " << fit_type << std::endl;
        outputFile << "Amp = " << Fit_Amp << " +/- " << Fit_Amp_e << std::endl;
        outputFile << "Mean = " << Fit_Mean << " +/- " << Fit_Mean_e << std::endl;
        outputFile << "Sigma = " << Fit_Sigma << " +/- " << Fit_Sigma_e << std::endl;
        outputFile << "---------------" << std::endl;
        outputFile.close(); // Close the file
    }else 
    {
        std::cerr << "Unable to open the file." << std::endl;
    }

    TLegend *legend = new TLegend(0.5, 0.6, 0.75, 0.9);
    legend->SetTextSize(0.025);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(2);
    legend->SetMargin(0.1); // Set a 10% margin

    std::stringstream amp_stream, amp_stream_error;
    amp_stream << std::fixed << std::setprecision(3) << Fit_Amp; // Set precision to 3 decimal places
    amp_stream_error << std::fixed << std::setprecision(3) << Fit_Amp_e; // Set precision to 3 decimal places
    std::string amp_legend = "Amp: " + amp_stream.str() + " +/- " + amp_stream_error.str();
    legend->AddEntry((TObject*)0, amp_legend.c_str(), "");

    std::stringstream mean_stream, mean_stream_error;
    mean_stream << std::fixed << std::setprecision(3) << Fit_Mean; // Set precision to 3 decimal places
    mean_stream_error << std::fixed << std::setprecision(3) << Fit_Mean_e; // Set precision to 3 decimal places
    std::string mean_legend = "Mean: " + mean_stream.str() + " +/- " + mean_stream_error.str();
    legend->AddEntry((TObject*)0, mean_legend.c_str(), "");

    std::stringstream sigma_stream, sigma_stream_error;
    sigma_stream << std::fixed << std::setprecision(3) << Fit_Sigma; // Set precision to 3 decimal places
    sigma_stream_error << std::fixed << std::setprecision(3) << Fit_Sigma_e; // Set precision to 3 decimal places
    std::string sigma_legend = "Sigma: " + sigma_stream.str() + " +/- " + sigma_stream_error.str();
    legend->AddEntry((TObject*)0, sigma_legend.c_str(), "");

    TCanvas* canvas = new TCanvas("canvas", "Fit Result Gaussian");
    if(log_style)
    {
        canvas->SetLogy(); 
    }
    hist->SetTitle(("Channel " + std::to_string(channel) + " with multiplicity " + std::to_string(MultiplicityCut)).c_str());
    hist->GetXaxis()->SetTitle("dt between pulse and trigger in ns");
    hist->GetYaxis()->SetTitle(("Counts/"+to_string(axis_label_value)+"ns").c_str());
    hist->Draw(); // Draw the original graph
    fit_function->Draw("same"); // Draw the fitted function on top
    legend->Draw();
    canvas->Update();
    std::cout << "Saving to " << m_data->TD.Path_Out << std::endl;
    std::string name = m_data->TD.Path_Out+"Gaussian_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas->SaveAs(name.c_str());  // Save the plot to a file
    canvas->Draw();
    
    // Clean up
    delete canvas;
    delete hist;
    delete fit_function;
}

void GetDecayTime::Use_Gaussian_plus_Exp2(int channel, TH1F *hist)
{
    float axis_label_value = (bins_end-bins_start)/n_bins;

    fit_function = new TF1("Gaussian", Gauss_plus_Exp2, fit_start, fit_end, 10);
    fit_function->SetParameters(guess_t0, guess_n0, guess_tau0, guess_n1, guess_tau1, guess_x0, guess_offset, guess_gauss_amp, guess_gauss_mean, guess_gauss_sigma);

    fit_function->SetParLimits(0,   0,      100);
    fit_function->SetParLimits(1,   0,      10*guess_n0);
    fit_function->SetParLimits(2,   0,      10*guess_tau0);
    fit_function->SetParLimits(3,   0,      10*guess_n1);
    fit_function->SetParLimits(4,   0,      10*guess_tau1);
    fit_function->SetParLimits(5,   -100,   100);
    fit_function->SetParLimits(6,   0,      10*guess_offset);

    fit_function->SetParLimits(7,   0,       10*guess_gauss_amp);
    fit_function->SetParLimits(8,   -100,    100);
    fit_function->SetParLimits(9,   0,       10);

    hist->Fit(fit_function, fit_options.c_str());

    Double_t Fit_t0 = fit_function->GetParameter(0);
    Double_t Fit_t0_e = fit_function->GetParError(0);

    Double_t Fit_n_0 = fit_function->GetParameter(1);
    Double_t Fit_n_0_e = fit_function->GetParError(1);

    Double_t Fit_tau_0 = fit_function->GetParameter(2);
    Double_t Fit_tau_0_e = fit_function->GetParError(2);

    Double_t Fit_n_1 = fit_function->GetParameter(3);
    Double_t Fit_n_1_e = fit_function->GetParError(3);

    Double_t Fit_tau_1 =  fit_function->GetParameter(4);
    Double_t Fit_tau_1_e = fit_function->GetParError(4);

    Double_t Fit_x0 = fit_function->GetParameter(5);
    Double_t Fit_x0_e = fit_function->GetParError(5);    

    Double_t Fit_Offset = fit_function->GetParameter(6);
    Double_t Fit_Offset_e = fit_function->GetParError(6);

    Double_t Fit_Amp = fit_function->GetParameter(7);
    Double_t Fit_Amp_e = fit_function->GetParError(7);

    Double_t Fit_Mean = fit_function->GetParameter(8);
    Double_t Fit_Mean_e = fit_function->GetParError(8);

    Double_t Fit_Sigma = fit_function->GetParameter(9);
    Double_t Fit_Sigma_e = fit_function->GetParError(9);

    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << std::endl;
    std::cout << "t0 = " << Fit_t0 << " +/- " << Fit_t0_e << " , x0 = " << Fit_x0 << " +/- " << Fit_x0_e << " , offset = " << Fit_Offset << " +/- " << Fit_Offset_e << std::endl;
    std::cout << BLUE << "Amp = " << Fit_Amp << " +/- " << Fit_Amp_e << RESET << std::endl;
    std::cout << BLUE << "Mean = " << Fit_Mean << " +/- " << Fit_Mean_e << RESET << std::endl;
    std::cout << GREEN << "Sigma = " << Fit_Sigma << " +/- " << Fit_Sigma_e << RESET << std::endl;
    printf("%s\n","---------------");
    std::cout << BLUE << "n0 = " << Fit_n_0 << " +/- " << Fit_n_0_e << RESET << ", " << GREEN << "tau0 = " << Fit_tau_0 << " +/- " << Fit_tau_0_e << RESET << std::endl;
    std::cout << BLUE << "n1 = " << Fit_n_1 << " +/- " << Fit_n_1_e << RESET << ", " << GREEN << "tau1 = " << Fit_tau_1 << " +/- " << Fit_tau_1_e << RESET << std::endl;
    printf("%s\n","---------------");

    std::string outputName = m_data->TD.Path_Out+"FitData_DecayTime.txt";
    std::ofstream outputFile(outputName.c_str(), std::ios_base::out | std::ios_base::app); // Open a file for writing
    if(outputFile.is_open()) 
    {
        // Write the formatted string to the file
        outputFile << "Fit result, channel " << channel << ", fit type was " << fit_type << std::endl;
        outputFile << "t0 = " << Fit_t0 << " +/- " << Fit_t0_e << " , x0 = " << Fit_x0 << " +/- " << Fit_x0_e << std::endl;
        outputFile << "Amp = " << Fit_Amp << " +/- " << Fit_Amp_e << std::endl;
        outputFile << "Mean = " << Fit_Mean << " +/- " << Fit_Mean_e << std::endl;
        outputFile << "Sigma = " << Fit_Sigma << " +/- " << Fit_Sigma_e << std::endl;
        outputFile << "---------------" << std::endl;
        outputFile << "n0 = " << Fit_n_0 << " +/- " << Fit_n_0_e << ", " << "tau0 = " << Fit_tau_0 << " +/- " << Fit_tau_0_e << std::endl;
        outputFile << "n1 = " << Fit_n_1 << " +/- " << Fit_n_1_e << ", " << "tau1 = " << Fit_tau_1 << " +/- " << Fit_tau_1_e << std::endl;
        outputFile << "---------------" << std::endl;
        outputFile.close(); // Close the file
    }else 
    {
        std::cerr << "Unable to open the file." << std::endl;
    }

    TLegend *legend = new TLegend(0.5, 0.6, 0.75, 0.9);
    legend->SetTextSize(0.025);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(2);
    legend->SetMargin(0.1); // Set a 10% margin

    std::stringstream t0_stream, t0_stream_error;
    t0_stream << std::fixed << std::setprecision(3) << Fit_t0; // Set precision to 3 decimal places
    t0_stream_error << std::fixed << std::setprecision(3) << Fit_t0_e; // Set precision to 3 decimal places
    std::string t0_legend = "t0: " + t0_stream.str() + " +/- " + t0_stream_error.str();
    legend->AddEntry((TObject*)0, t0_legend.c_str(), "");

    std::stringstream n0_stream, n0_stream_error;
    n0_stream << std::fixed << std::setprecision(3) << Fit_n_0; // Set precision to 3 decimal places
    n0_stream_error << std::fixed << std::setprecision(3) << Fit_n_0_e; // Set precision to 3 decimal places
    std::string n0_legend = "n0: " + n0_stream.str() + " +/- " + n0_stream_error.str();
    legend->AddEntry((TObject*)0, n0_legend.c_str(), "");

    std::stringstream tau0_stream, tau0_stream_error;
    tau0_stream << std::fixed << std::setprecision(3) << Fit_tau_0; // Set precision to 3 decimal places
    tau0_stream_error << std::fixed << std::setprecision(3) << Fit_tau_0_e; // Set precision to 3 decimal places
    std::string tau0_legend = "tau0: " + tau0_stream.str() + " +/- " + tau0_stream_error.str();
    legend->AddEntry((TObject*)0, tau0_legend.c_str(), "");

    std::stringstream n1_stream, n1_stream_error;
    n1_stream << std::fixed << std::setprecision(3) << Fit_n_1; // Set precision to 3 decimal places
    n1_stream_error << std::fixed << std::setprecision(3) << Fit_n_1_e; // Set precision to 3 decimal places
    std::string n1_legend = "n1: " + n1_stream.str() + " +/- " + n1_stream_error.str();
    legend->AddEntry((TObject*)0, n1_legend.c_str(), "");

    std::stringstream tau1_stream, tau1_stream_error;
    tau1_stream << std::fixed << std::setprecision(3) << Fit_tau_1; // Set precision to 3 decimal places
    tau1_stream_error << std::fixed << std::setprecision(3) << Fit_tau_1_e; // Set precision to 3 decimal places
    std::string tau1_legend = "tau1: " + tau1_stream.str() + " +/- " + tau1_stream_error.str();
    legend->AddEntry((TObject*)0, tau1_legend.c_str(), "");

    std::stringstream amp_stream, amp_stream_error;
    amp_stream << std::fixed << std::setprecision(3) << Fit_Amp; // Set precision to 3 decimal places
    amp_stream_error << std::fixed << std::setprecision(3) << Fit_Amp_e; // Set precision to 3 decimal places
    std::string amp_legend = "Amp: " + amp_stream.str() + " +/- " + amp_stream_error.str();
    legend->AddEntry((TObject*)0, amp_legend.c_str(), "");

    std::stringstream mean_stream, mean_stream_error;
    mean_stream << std::fixed << std::setprecision(3) << Fit_Mean; // Set precision to 3 decimal places
    mean_stream_error << std::fixed << std::setprecision(3) << Fit_Mean_e; // Set precision to 3 decimal places
    std::string mean_legend = "Mean: " + mean_stream.str() + " +/- " + mean_stream_error.str();
    legend->AddEntry((TObject*)0, mean_legend.c_str(), "");

    std::stringstream sigma_stream, sigma_stream_error;
    sigma_stream << std::fixed << std::setprecision(3) << Fit_Sigma; // Set precision to 3 decimal places
    sigma_stream_error << std::fixed << std::setprecision(3) << Fit_Sigma_e; // Set precision to 3 decimal places
    std::string sigma_legend = "Sigma: " + sigma_stream.str() + " +/- " + sigma_stream_error.str();
    legend->AddEntry((TObject*)0, sigma_legend.c_str(), "");

    std::stringstream label;
    label << std::fixed << std::setprecision(2) << axis_label_value; // Set precision to 3 decimal places

    singleExpFunc = new TF1("singleExpFunc", Sum_of_Exp_Normal_2, bins_start, bins_end, 6); 
    singleGaussFunc = new TF1("singleGaussFunc", Gauss1, bins_start, bins_end, 3); 
    linearFunc0 = new TF1("linearFunc0", LinearF0, bins_start, bins_end, 2); 

    singleExpFunc->SetParameters(Fit_t0, Fit_n_0, Fit_tau_0, Fit_n_1, Fit_tau_1, Fit_x0);
    singleGaussFunc->SetParameters(Fit_Amp, Fit_Mean, Fit_Sigma);
    linearFunc0->SetParameter(0,Fit_Offset);

    singleExpFunc->SetLineColor(kViolet);
    singleExpFunc->SetLineStyle(2);
    singleGaussFunc->SetLineColor(kGreen);
    singleGaussFunc->SetLineStyle(2);
    linearFunc0->SetLineColor(kYellow);
    linearFunc0->SetLineStyle(2);

    TCanvas* canvas = new TCanvas("canvas", "Fit Result Gaussian");
    if(log_style)
    {
        canvas->SetLogy(); 
    }
    hist->SetTitle(("Channel " + std::to_string(channel) + " with multiplicity " + std::to_string(MultiplicityCut)).c_str());
    hist->GetXaxis()->SetTitle("dt between pulse and trigger in ns");
    hist->GetYaxis()->SetTitle((std::string("Counts/")+label.str()+"ns").c_str());
    hist->Draw(); // Draw the original graph
    fit_function->Draw("same"); // Draw the fitted function on top
    singleExpFunc->Draw("same");
    singleGaussFunc->Draw("same");
    linearFunc0->Draw("same");
    legend->Draw();
    canvas->Update();
    std::string name = m_data->TD.Path_Out+"Gaussian_p_Exp2_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas->SaveAs(name.c_str());  // Save the plot to a file
    canvas->Draw();
    
    // Clean up
    delete canvas;
    delete hist;
    delete fit_function;
}

void GetDecayTime::InitRoot()
{
    m_data->TD.TTree_Analysis_DecayTime = new TTree("TTree_Analysis_DecayTime", "TTree_Analysis_DecayTime");
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch0", &Decay_Ch0);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch1", &Decay_Ch1);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch2", &Decay_Ch2);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch3", &Decay_Ch3);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch4", &Decay_Ch4);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch5", &Decay_Ch5);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch6", &Decay_Ch6);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch7", &Decay_Ch7);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch8", &Decay_Ch8);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch9", &Decay_Ch9);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch10", &Decay_Ch10);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch11", &Decay_Ch11);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch12", &Decay_Ch12);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch13", &Decay_Ch13);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch14", &Decay_Ch14);
    m_data->TD.TTree_Analysis_DecayTime->Branch("Decay_Ch15", &Decay_Ch15);  
}



void GetDecayTime::Use_Fold4(int channel, TH1F *hist)
{
    float axis_label_value = (bins_end-bins_start)/n_bins;

    cout << "Before Fit" << endl;
    fit_function = new TF1("fit_function", Sum_of_Exp_Fold4, fit_start, fit_end, 9); 
    cout << "After Fit" << endl;
    fit_function->SetParameters(guess_t0, guess_n0, guess_tau0, guess_n1, guess_tau1, guess_n2, guess_tau2, guess_n3, guess_tau3);
    m_data->FitModel.summation_order = 4;
    
    // Set bounds for x0 between -100 and 100
    fit_function->SetParLimits(0, 1.5*guess_t0, 0.5*guess_t0);
    // Set bounds for n_0 between 0 and 100
    fit_function->SetParLimits(1, 0, 10*guess_n0);
    // Set bounds for tau_0 between 0 and 10
    fit_function->SetParLimits(2, 0, 10*guess_tau0);
    // Set bounds for n_1 between 0 and 100
    fit_function->SetParLimits(3, 0, 10*guess_n1);
    // Set bounds for tau_1 between 0 and 10
    fit_function->SetParLimits(4, 0, 10*guess_tau1);
    // Set bounds for n_2 between 0 and 100
    fit_function->SetParLimits(5, 0, 10*guess_n2);
    // Set bounds for tau_2 between 0 and 10
    fit_function->SetParLimits(6, 0, 10*guess_tau2);
    // Set bounds for n_3 between 0 and 100
    fit_function->SetParLimits(7, 0, 10*guess_n2);
    // Set bounds for tau_3 between 0 and 10
    fit_function->SetParLimits(8, 0, 10*guess_tau2);

    // Perform the fit
    hist->Fit(fit_function, fit_options.c_str()); // "S" for using the sum of squares method

    Double_t fit_t0 = fit_function->GetParameter(0); 
    Double_t fit_t0_error = fit_function->GetParError(0); 

    Double_t fit_n_0 = fit_function->GetParameter(1);
    Double_t fit_n_0_error = fit_function->GetParError(1);

    Double_t fit_tau_0 = fit_function->GetParameter(2);
    Double_t fit_tau_0_error = fit_function->GetParError(2);

    Double_t fit_n_1 = fit_function->GetParameter(3);
    Double_t fit_n_1_error = fit_function->GetParError(3);

    Double_t fit_tau_1 =  fit_function->GetParameter(4);
    Double_t fit_tau_1_error = fit_function->GetParError(4);

    Double_t fit_n_2 = fit_function->GetParameter(5);
    Double_t fit_n_2_error = fit_function->GetParError(5);

    Double_t fit_tau_2 = fit_function->GetParameter(6);
    Double_t fit_tau_2_error = fit_function->GetParError(6);

    Double_t fit_n_3 = fit_function->GetParameter(7);
    Double_t fit_n_3_error = fit_function->GetParError(7);

    Double_t fit_tau_3 = fit_function->GetParameter(8);
    Double_t fit_tau_3_error = fit_function->GetParError(8);
    

    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << ": t0 = " << fit_t0 << " +/- " << fit_t0_error <<  std::endl;
    std::cout << BLUE << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << RESET << ", " << GREEN << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << RESET << std::endl;
    std::cout << BLUE << "n1 = " << fit_n_1 << " +/- " << fit_n_1_error << RESET << ", " << GREEN << "tau1 = " << fit_tau_1 << " +/- " << fit_tau_1_error << RESET << std::endl;
    std::cout << BLUE << "n2 = " << fit_n_2 << " +/- " << fit_n_2_error << RESET << ", " << GREEN << "tau2 = " << fit_tau_2 << " +/- " << fit_tau_2_error << RESET << std::endl;
    std::cout << BLUE << "n3 = " << fit_n_3 << " +/- " << fit_n_3_error << RESET << ", " << GREEN << "tau3 = " << fit_tau_3 << " +/- " << fit_tau_3_error << RESET << std::endl;
    printf("%s\n","---------------");


    // Create a legend
    TLegend *legend = new TLegend(0.5, 0.6, 0.75, 0.9);
    legend->SetTextSize(0.025);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(2);
    legend->SetMargin(0.1); // Set a 10% margin

    // Add fit parameters and errors to the legend
    AddLegend(legend, fit_t0, fit_t0_error, "t0");
    AddLegend(legend, fit_n_0, fit_n_0_error, "n0");
    AddLegend(legend, fit_tau_0, fit_tau_0_error, "tau0");
    AddLegend(legend, fit_n_1, fit_n_1_error, "n1");
    AddLegend(legend, fit_tau_1, fit_tau_1_error, "tau1");
    AddLegend(legend, fit_n_2, fit_n_2_error, "n2");
    AddLegend(legend, fit_tau_2, fit_tau_2_error, "tau2");
    AddLegend(legend, fit_n_3, fit_n_3_error, "n3");
    AddLegend(legend, fit_tau_3, fit_tau_3_error, "tau3");

    std::stringstream label;
    label << std::fixed << std::setprecision(2) << axis_label_value; // Set precision to 3 decimal places

    singleExpFunc0 = new TF1("singleExpFunc0", [this](double *x, double *par)->Double_t{return m_data->FitModel.ExponentialDecay_ConvSingle(x, par);}, bins_start, bins_end, 3); // 4 parameters for Single_Exp
    singleExpFunc1 = new TF1("singleExpFunc1", [this](double *x, double *par)->Double_t{return m_data->FitModel.ExponentialDecay_ConvSingle(x, par);}, bins_start, bins_end, 3); // 4 parameters for Single_Exp
    singleExpFunc2 = new TF1("singleExpFunc2", [this](double *x, double *par)->Double_t{return m_data->FitModel.ExponentialDecay_ConvSingle(x, par);}, bins_start, bins_end, 3); // 4 parameters for Single_Exp
    singleExpFunc3 = new TF1("singleExpFunc3", [this](double *x, double *par)->Double_t{return m_data->FitModel.ExponentialDecay_ConvSingle(x, par);}, bins_start, bins_end, 3); // 4 parameters for Single_Exp

    singleExpFunc0->SetParameters(fit_n_0, fit_tau_0, fit_t0);
    singleExpFunc1->SetParameters(fit_n_1, fit_tau_1, fit_t0);
    singleExpFunc2->SetParameters(fit_n_2, fit_tau_2, fit_t0);
    singleExpFunc3->SetParameters(fit_n_3, fit_tau_3, fit_t0);

    singleExpFunc0->SetLineColor(kGray);
    singleExpFunc0->SetLineStyle(2);
    singleExpFunc1->SetLineColor(kGreen);
    singleExpFunc1->SetLineStyle(2);
    singleExpFunc2->SetLineColor(kViolet);
    singleExpFunc2->SetLineStyle(2);
    singleExpFunc3->SetLineColor(kYellow);
    singleExpFunc3->SetLineStyle(2);

    TCanvas* canvas = new TCanvas("canvas", "Fit Result", 800, 600);
    if(log_style)
    {
        canvas->SetLogy(); 
    }
    hist->SetTitle(("Channel " + std::to_string(channel) + " with multiplicity " + std::to_string(MultiplicityCut)).c_str());
    hist->GetXaxis()->SetTitle("dt between pulse and trigger in ns");
    hist->GetYaxis()->SetTitle((std::string("Counts/")+label.str()+"ns").c_str());
    hist->Draw(); // Draw the original graph
    fit_function->Draw("same"); // Draw the fitted function on top
    singleExpFunc0->Draw("same");
    singleExpFunc1->Draw("same");
    singleExpFunc2->Draw("same");
    singleExpFunc3->Draw("same");
    legend->Draw();
    canvas->Update();
    std::string name = m_data->TD.Path_Out+"Fold4_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas->SaveAs(name.c_str());  // Save the plot to a file
    canvas->Draw();
    delete canvas;

    // Clean up
    delete hist;
    delete fit_function;
}


    // // Set bounds for t0 between -100 and 100
    // fit_function->SetParLimits(0, 1.5*guess_t0, 0.5*guess_t0);
    // // Set bounds for tauR between -100 and 100
    // fit_function->SetParLimits(1, 0, 10*guess_tauR);
    // // Set bounds for n_0 between 0 and 100
    // fit_function->SetParLimits(2, 0, 10*guess_n0);
    // // Set bounds for tau_0 between 0 and 10
    // fit_function->SetParLimits(3, 0, 10*guess_tau0);
    // // Set bounds for n_1 between 0 and 100
    // fit_function->SetParLimits(4, 0, 10*guess_n1);
    // // Set bounds for tau_1 between 0 and 10
    // fit_function->SetParLimits(5, 0, 10*guess_tau1);
    // // Set bounds for n_2 between 0 and 100
    // fit_function->SetParLimits(6, 0, 10*guess_n2);
    // // Set bounds for tau_2 between 0 and 10
    // fit_function->SetParLimits(7, 0, 10*guess_tau2);
    // // Set bounds for n_3 between 0 and 100
    // fit_function->SetParLimits(8, 0, 10*guess_n3);
    // // Set bounds for tau_3 between 0 and 10
    // fit_function->SetParLimits(9, 0, 10*guess_tau3);

    // Perform the fit
    
    //hist->Fit("Convolution", fit_options.c_str()); // "S" for using the sum of squares method
/*
    Double_t fit_t0 = fit_function->GetParameter(0);
    Double_t fit_t0_error = fit_function->GetParError(0);

    Double_t fit_tauR = fit_function->GetParameter(1);
    Double_t fit_tauR_error = fit_function->GetParError(1); 

    Double_t fit_n_0 = fit_function->GetParameter(2);
    Double_t fit_n_0_error = fit_function->GetParError(2);

    Double_t fit_tau_0 = fit_function->GetParameter(3);
    Double_t fit_tau_0_error = fit_function->GetParError(3);

    Double_t fit_n_1 = fit_function->GetParameter(4);
    Double_t fit_n_1_error = fit_function->GetParError(4);

    Double_t fit_tau_1 =  fit_function->GetParameter(5);
    Double_t fit_tau_1_error = fit_function->GetParError(5);

    Double_t fit_n_2 = fit_function->GetParameter(6);
    Double_t fit_n_2_error = fit_function->GetParError(6);

    Double_t fit_tau_2 = fit_function->GetParameter(7);
    Double_t fit_tau_2_error = fit_function->GetParError(7);

    Double_t fit_n_3 = fit_function->GetParameter(8);
    Double_t fit_n_3_error = fit_function->GetParError(8);

    Double_t fit_tau_3 = fit_function->GetParameter(9);
    Double_t fit_tau_3_error = fit_function->GetParError(9);
    

    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << ": t0 = " << fit_t0 << " +/- " << fit_t0_error << " , tauR = " << fit_tauR << " +/- " << fit_tauR_error << std::endl;
    std::cout << BLUE << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << RESET << ", " << GREEN << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << RESET << std::endl;
    std::cout << BLUE << "n1 = " << fit_n_1 << " +/- " << fit_n_1_error << RESET << ", " << GREEN << "tau1 = " << fit_tau_1 << " +/- " << fit_tau_1_error << RESET << std::endl;
    std::cout << BLUE << "n2 = " << fit_n_2 << " +/- " << fit_n_2_error << RESET << ", " << GREEN << "tau2 = " << fit_tau_2 << " +/- " << fit_tau_2_error << RESET << std::endl;
    std::cout << BLUE << "n3 = " << fit_n_3 << " +/- " << fit_n_3_error << RESET << ", " << GREEN << "tau3 = " << fit_tau_3 << " +/- " << fit_tau_3_error << RESET << std::endl;
    printf("%s\n","---------------");

    // Create a legend
    TLegend *legend = new TLegend(0.5, 0.6, 0.75, 0.9);
    legend->SetTextSize(0.025);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(2);
    legend->SetMargin(0.1); // Set a 10% margin

    // Add fit parameters and errors to the legend

    AddLegend(legend, fit_t0, fit_t0_error, "t0");
    AddLegend(legend, fit_tauR, fit_tauR_error, "tauR");
    AddLegend(legend, fit_n_0, fit_n_0_error, "n0");
    AddLegend(legend, fit_tau_0, fit_tau_0_error, "tau0");
    AddLegend(legend, fit_n_1, fit_n_1_error, "n1");
    AddLegend(legend, fit_tau_1, fit_tau_1_error, "tau1");
    AddLegend(legend, fit_n_2, fit_n_2_error, "n2");
    AddLegend(legend, fit_tau_2, fit_tau_2_error, "tau2");
    AddLegend(legend, fit_n_3, fit_n_3_error, "n3");
    AddLegend(legend, fit_tau_3, fit_tau_3_error, "tau3");

    std::stringstream label;
    label << std::fixed << std::setprecision(2) << axis_label_value; // Set precision to 3 decimal places

    singleExpFunc0 = new TF1("singleExpFunc0", [this](double *x, double *par)->Double_t{return m_data->FitModel.ExponentialDecay_RiseSum(x, par);}, fit_start, fit_end, 4); // 4 parameters for Single_Exp
    singleExpFunc1 = new TF1("singleExpFunc1", [this](double *x, double *par)->Double_t{return m_data->FitModel.ExponentialDecay_RiseSum(x, par);}, fit_start, fit_end, 4); // 4 parameters for Single_Exp
    singleExpFunc2 = new TF1("singleExpFunc2", [this](double *x, double *par)->Double_t{return m_data->FitModel.ExponentialDecay_RiseSum(x, par);}, fit_start, fit_end, 4); // 4 parameters for Single_Exp
    singleExpFunc3 = new TF1("singleExpFunc3", [this](double *x, double *par)->Double_t{return m_data->FitModel.ExponentialDecay_RiseSum(x, par);}, fit_start, fit_end, 4); // 4 parameters for Single_Exp

    singleExpFunc0->SetParameters(fit_t0, fit_tauR, fit_n_0, fit_tau_0);
    singleExpFunc1->SetParameters(fit_t0, fit_tauR, fit_n_1, fit_tau_1);
    singleExpFunc2->SetParameters(fit_t0, fit_tauR, fit_n_2, fit_tau_2);
    singleExpFunc3->SetParameters(fit_t0, fit_tauR, fit_n_3, fit_tau_3);   

    singleExpFunc0->SetLineColor(kBlue);
    singleExpFunc0->SetLineStyle(2);
    singleExpFunc1->SetLineColor(kGreen);
    singleExpFunc1->SetLineStyle(2);
    singleExpFunc2->SetLineColor(kViolet);
    singleExpFunc2->SetLineStyle(2);
    singleExpFunc3->SetLineColor(kGray);
    singleExpFunc3->SetLineStyle(2);

    TCanvas* canvas = new TCanvas("canvas", "Fit Result");
    if(log_style)
    {
        canvas->SetLogy(); 
    }
    hist->SetTitle(("Channel " + std::to_string(channel) + " with multiplicity " + std::to_string(MultiplicityCut)).c_str());
    hist->GetXaxis()->SetTitle("dt between pulse and trigger in ns");
    hist->GetYaxis()->SetTitle((std::string("Counts/")+label.str()+"ns").c_str());
    hist->Draw(); // Draw the original graph
    fit_function->Draw("same"); // Draw the fitted function on top
    singleExpFunc0->Draw("same");
    singleExpFunc1->Draw("same");
    singleExpFunc2->Draw("same");
    singleExpFunc3->Draw("same");
    legend->Draw();
    canvas->Update();
    std::string name = m_data->TD.Path_Out+"Exp4_ch_"+to_string(channel)+"_C"+to_string(MultiplicityCut)+".png";
    canvas->SaveAs(name.c_str());  // Save the plot to a file
    canvas->Draw();
    delete canvas;

    // Clean up
    delete hist;
    delete fit_function;*/