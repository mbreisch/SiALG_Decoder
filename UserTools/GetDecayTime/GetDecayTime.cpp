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

    if(!m_variables.Get("gu",threshold_multiplier)) threshold_multiplier=0.5;      
    if(!m_variables.Get("guess_x0",guess_x0)) guess_x0=0.5;     
    if(!m_variables.Get("guess_t0",guess_t0)) guess_t0=0.5;     
    if(!m_variables.Get("guess_n0",guess_n0)) guess_n0=0.5;     
    if(!m_variables.Get("guess_tau0",guess_tau0)) guess_tau0=0.5;     
    if(!m_variables.Get("guess_n1",guess_n1)) guess_n1=0.5;     
    if(!m_variables.Get("guess_tau1",guess_tau1)) guess_tau1=0.5;     
    if(!m_variables.Get("guess_n2",guess_n2)) guess_n2=0.5;     
    if(!m_variables.Get("guess_tau2",guess_tau2)) guess_tau2=0.5;    
    if(!m_variables.Get("guess_n3",guess_n3)) guess_n3=0.5;     
    if(!m_variables.Get("guess_tau3",guess_tau3)) guess_tau3=0.5;  

    if(!m_variables.Get("guess_gauss_amp",guess_gauss_amp)) guess_gauss_amp=0.5;    
    if(!m_variables.Get("guess_gauss_mean",guess_gauss_mean)) guess_gauss_mean=0.5;     
    if(!m_variables.Get("guess_gauss_sigma",guess_gauss_sigma)) guess_gauss_sigma=0.5;      

    if(!m_variables.Get("fit_type",fit_type)) fit_type="normal";    
    if(!m_variables.Get("fit_options",fit_options)) fit_options="RQ";    
    if(!m_variables.Get("log_style",log_style)) log_style=true;    

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

Double_t Sum_of_Exp_Normal_4(Double_t *x, Double_t *par) 
{
    Double_t tau_R = par[0];
    Double_t n_0 = par[1];
    Double_t tau_0 = par[2];
    Double_t n_1 = par[3];
    Double_t tau_1 = par[4];
    Double_t n_2 = par[5];
    Double_t tau_2 = par[6];
    Double_t n_3 = par[7];
    Double_t tau_3 = par[8];
    Double_t x_0 = par[9];

    return n_0 * (  (   TMath::Exp(-(x[0]-x_0)/tau_0) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_0 - tau_R)   ) +
            n_1 * (  (   TMath::Exp(-(x[0]-x_0)/tau_1) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_1 - tau_R)   ) +
            n_2 * (  (   TMath::Exp(-(x[0]-x_0)/tau_2) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_2 - tau_R)   ) +
            n_3 * (  (   TMath::Exp(-(x[0]-x_0)/tau_3) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_3 - tau_R)   );
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

    return n_0 * (  (   TMath::Exp(-(x[0]-x_0)/tau_0) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_0 - tau_R)   ) +
            n_1 * (  (   TMath::Exp(-(x[0]-x_0)/tau_1) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_1 - tau_R)   ) +
            n_2 * (  (   TMath::Exp(-(x[0]-x_0)/tau_2) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_2 - tau_R)   );
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
    Double_t tau_R = par[0];
    Double_t n_0 = par[1];
    Double_t tau_0 = par[2];
    Double_t x_0 = par[3];

    return n_0 * (  (   TMath::Exp(-(x[0]-x_0)/tau_0) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_0 - tau_R)   );
}

Double_t Gauss1(Double_t *x, Double_t *par)
{
    Double_t amp = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];

    return (amp / (TMath::Sqrt(2.0 * TMath::Pi())*sigma)) * TMath::Exp(-0.5 * TMath::Power((x[0] - mean) / sigma, 2));
}

Double_t Gauss_plus_Exp2(Double_t *x, Double_t *par)
{
    Double_t tau_R = par[0];
    Double_t n_0 = par[1];
    Double_t tau_0 = par[2];
    Double_t n_1 = par[3];
    Double_t tau_1 = par[4];
    Double_t x_0 = par[5];

    Double_t amp = par[6];
    Double_t mean = par[7];
    Double_t sigma = par[8];

    return (amp / (TMath::Sqrt(2.0 * TMath::Pi())*sigma)) * TMath::Exp(-0.5 * TMath::Power((x[0] - mean) / sigma, 2)) +
                (n_0 * (  (   TMath::Exp(-(x[0]-x_0)/tau_0) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_0 - tau_R)   ) +
                 n_1 * (  (   TMath::Exp(-(x[0]-x_0)/tau_1) - TMath::Exp(-(x[0]-x_0)/tau_R) ) /   (tau_1 - tau_R)   ));
}

void GetDecayTime::FitDecay(int channel, vector<float> data)
{
    if(data.size()==0){return;}

    TH1F *hist = new TH1F("hist", "Example Histogram", n_bins,bins_start,bins_end);
    for(float k: data) 
    {
        hist->Fill(k);
    }
    hist->Scale(1.0 / hist->Integral());

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
    }else 
    {
        std::cout << RED << "No Fit Type was chosen" << RESET << std::endl;
    }

}

void GetDecayTime::Use_Exp1(int channel, TH1F *hist)
{
    float axis_label_value = (bins_end-bins_start)/n_bins;

    fit_function = new TF1("fit_function", Single_Exp_Normal, fit_start, fit_end, 4); 
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

    fit_function = new TF1("fit_function", Sum_of_Exp_Normal_3, fit_start, fit_end, 8); 
    fit_function->SetParameters(guess_t0, guess_n0, guess_tau0, guess_n1, guess_tau1, guess_n2, guess_tau2, guess_x0);

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

    singleExpFunc0 = new TF1("singleExpFunc0", Single_Exp_Normal, fit_start, fit_end, 4); // 4 parameters for Single_Exp
    singleExpFunc1 = new TF1("singleExpFunc1", Single_Exp_Normal, fit_start, fit_end, 4); // 4 parameters for Single_Exp
    singleExpFunc2 = new TF1("singleExpFunc2", Single_Exp_Normal, fit_start, fit_end, 4); // 4 parameters for Single_Exp

    singleExpFunc0->SetParameters(fit_t0, fit_n_0, fit_tau_0, fit_x0);
    singleExpFunc1->SetParameters(fit_t0, fit_n_1, fit_tau_1, fit_x0);
    singleExpFunc2->SetParameters(fit_t0, fit_n_2, fit_tau_2, fit_x0);

    singleExpFunc0->SetLineColor(kGray);
    singleExpFunc0->SetLineStyle(2);
    singleExpFunc1->SetLineColor(kGreen);
    singleExpFunc1->SetLineStyle(2);
    singleExpFunc2->SetLineColor(kViolet);
    singleExpFunc2->SetLineStyle(2);

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

    fit_function = new TF1("fit_function", Sum_of_Exp_Normal_4, fit_start, fit_end, 10); // 8 is the number of parameters
    fit_function->SetParameters(guess_t0, guess_n0, guess_tau0, guess_n1, guess_tau1, guess_n2, guess_tau2, guess_n3, guess_tau3, guess_x0);

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
    // Set bounds for n_3 between 0 and 100
    fit_function->SetParLimits(7, 0, 10*guess_n3);
    // Set bounds for tau_3 between 0 and 10
    fit_function->SetParLimits(8, 0, 10*guess_tau3);
    // Set bounds for x_0 between 0 and 10
    fit_function->SetParLimits(9, -100, 100);

    // Perform the fit
    hist->Fit(fit_function, fit_options.c_str()); // "S" for using the sum of squares method

    Double_t fit_t0 = fit_function->GetParameter(0);
    Double_t fit_n_0 = fit_function->GetParameter(1);
    Double_t fit_tau_0 = fit_function->GetParameter(2);
    Double_t fit_n_1 = fit_function->GetParameter(3);
    Double_t fit_tau_1 =  fit_function->GetParameter(4);
    Double_t fit_n_2 = fit_function->GetParameter(5);
    Double_t fit_tau_2 = fit_function->GetParameter(6);
    Double_t fit_n_3 = fit_function->GetParameter(7);
    Double_t fit_tau_3 = fit_function->GetParameter(8);
    Double_t fit_x0 = fit_function->GetParameter(9); 

    Double_t fit_t0_error = fit_function->GetParError(0);
    Double_t fit_n_0_error = fit_function->GetParError(1);
    Double_t fit_tau_0_error = fit_function->GetParError(2);
    Double_t fit_n_1_error = fit_function->GetParError(3);
    Double_t fit_tau_1_error = fit_function->GetParError(4);
    Double_t fit_n_2_error = fit_function->GetParError(5);
    Double_t fit_tau_2_error = fit_function->GetParError(6);
    Double_t fit_n_3_error = fit_function->GetParError(7);
    Double_t fit_tau_3_error = fit_function->GetParError(8);
    Double_t fit_x0_error = fit_function->GetParError(9); 

    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << ": t0 = " << fit_t0 << " +/- " << fit_t0_error << " , x0 = " << fit_x0 << " +/- " << fit_x0_error << std::endl;
    std::cout << BLUE << "n0 = " << fit_n_0 << " +/- " << fit_n_0_error << RESET << ", " << GREEN << "tau0 = " << fit_tau_0 << " +/- " << fit_tau_0_error << RESET << std::endl;
    std::cout << BLUE << "n1 = " << fit_n_1 << " +/- " << fit_n_1_error << RESET << ", " << GREEN << "tau1 = " << fit_tau_1 << " +/- " << fit_tau_1_error << RESET << std::endl;
    std::cout << BLUE << "n2 = " << fit_n_2 << " +/- " << fit_n_2_error << RESET << ", " << GREEN << "tau2 = " << fit_tau_2 << " +/- " << fit_tau_2_error << RESET << std::endl;
    std::cout << BLUE << "n3 = " << fit_n_3 << " +/- " << fit_n_3_error << RESET << ", " << GREEN << "tau3 = " << fit_tau_3 << " +/- " << fit_tau_3_error << RESET << std::endl;
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
        outputFile << "n3 = " << fit_n_3 << " +/- " << fit_n_3_error << ", " << "tau3 = " << fit_tau_3 << " +/- " << fit_tau_3_error << std::endl;
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

    std::stringstream n3_stream, n3_stream_error;
    n3_stream << std::fixed << std::setprecision(3) << fit_n_3; // Set precision to 3 decimal places
    n3_stream_error << std::fixed << std::setprecision(3) << fit_n_3_error; // Set precision to 3 decimal places
    std::string n3_legend = "n3: " + n3_stream.str() + " +/- " + n3_stream_error.str();
    legend->AddEntry((TObject*)0, n3_legend.c_str(), "");

    std::stringstream tau3_stream, tau3_stream_error;
    tau3_stream << std::fixed << std::setprecision(3) << fit_tau_3; // Set precision to 3 decimal places
    tau3_stream_error << std::fixed << std::setprecision(3) << fit_tau_3_error; // Set precision to 3 decimal places
    std::string tau3_legend = "tau3: " + tau3_stream.str() + " +/- " + tau3_stream_error.str();
    legend->AddEntry((TObject*)0, tau3_legend.c_str(), "");

    std::stringstream label;
    label << std::fixed << std::setprecision(2) << axis_label_value; // Set precision to 3 decimal places

    singleExpFunc0 = new TF1("singleExpFunc0", Single_Exp_Normal, fit_start, fit_end, 4); // 4 parameters for Single_Exp
    singleExpFunc1 = new TF1("singleExpFunc1", Single_Exp_Normal, fit_start, fit_end, 4); // 4 parameters for Single_Exp
    singleExpFunc2 = new TF1("singleExpFunc2", Single_Exp_Normal, fit_start, fit_end, 4); // 4 parameters for Single_Exp
    singleExpFunc3 = new TF1("singleExpFunc3", Single_Exp_Normal, fit_start, fit_end, 4); // 4 parameters for Single_Exp

    singleExpFunc0->SetParameters(fit_t0, fit_n_0, fit_tau_0, fit_x0);
    singleExpFunc1->SetParameters(fit_t0, fit_n_1, fit_tau_1, fit_x0);
    singleExpFunc2->SetParameters(fit_t0, fit_n_2, fit_tau_2, fit_x0);
    singleExpFunc3->SetParameters(fit_t0, fit_n_3, fit_tau_3, fit_x0);   

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
    delete fit_function;
}

void GetDecayTime::Use_Gaussian(int channel, TH1F *hist)
{
    int axis_label_value = (bins_end-bins_start)/n_bins;

    fit_function = new TF1("Gaussian", Gauss1, fit_start, fit_end, 3);
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

    fit_function = new TF1("Gaussian", Gauss_plus_Exp2, fit_start, fit_end, 9);
    fit_function->SetParameters(guess_t0, guess_n0, guess_tau0, guess_n1, guess_tau1, guess_x0, guess_gauss_amp, guess_gauss_mean, guess_gauss_sigma);

    fit_function->SetParLimits(0,   0,      100);
    fit_function->SetParLimits(1,   0,      1);
    fit_function->SetParLimits(2,   0,      10*guess_tau0);
    fit_function->SetParLimits(3,   0,      1);
    fit_function->SetParLimits(4,   0,      10*guess_tau1);
    fit_function->SetParLimits(5,   -100,   100);

    fit_function->SetParLimits(6,   0,       10*guess_gauss_amp);
    fit_function->SetParLimits(7,   -100,    100);
    fit_function->SetParLimits(8,   0,       10);

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

    Double_t Fit_Amp = fit_function->GetParameter(6);
    Double_t Fit_Amp_e = fit_function->GetParError(6);

    Double_t Fit_Mean = fit_function->GetParameter(7);
    Double_t Fit_Mean_e = fit_function->GetParError(7);

    Double_t Fit_Sigma = fit_function->GetParameter(8);
    Double_t Fit_Sigma_e = fit_function->GetParError(8);

    // Print the results
    printf("%s\n","---------------");
    std::cout << "Fit result, channel " << channel << std::endl;
    std::cout << "t0 = " << Fit_t0 << " +/- " << Fit_t0_e << " , x0 = " << Fit_x0 << " +/- " << Fit_x0_e << std::endl;
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

    singleExpFunc = new TF1("singleExpFunc", Sum_of_Exp_Normal_2, fit_start, fit_end, 6); 
    singleGaussFunc = new TF1("singleGaussFunc", Gauss1, fit_start, fit_end, 3); 

    singleExpFunc->SetParameters(Fit_t0, Fit_n_0, Fit_tau_0, Fit_n_1, Fit_tau_1, Fit_x0);
    singleGaussFunc->SetParameters(Fit_Amp, Fit_Mean, Fit_Sigma);

    singleExpFunc->SetLineColor(kViolet);
    singleExpFunc->SetLineStyle(2);
    singleGaussFunc->SetLineColor(kGreen);
    singleGaussFunc->SetLineStyle(2);

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
    Double_t n_3 = par[8];
    Double_t tau_3 = par[9];

    return TMath::Exp(-(x[0]*x[0])/(sigma*sigma)) * 
            (
                (n_0/tau_0) * TMath::Exp(-(x[0] - x0)/tau_0) +
                (n_1/tau_1) * TMath::Exp(-(x[0] - x0)/tau_1) +
                (n_2/tau_2) * TMath::Exp(-(x[0] - x0)/tau_2) +
                (n_3/tau_3) * TMath::Exp(-(x[0] - x0)/tau_3)
            );
}
Double_t Single_Exp_Fold(Double_t *x, Double_t *par)
{
    Double_t sigma = par[0];
    Double_t x0 = par[1];
    Double_t n_0 = par[2];
    Double_t tau_0 = par[3];

    return  TMath::Exp(-(x[0]*x[0])/(sigma*sigma))* 
            (
                (n_0/tau_0) * TMath::Exp(-(x[0] - x0)/tau_0)
            );
}