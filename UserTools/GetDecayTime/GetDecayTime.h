#ifndef GetDecayTime_H
#define GetDecayTime_H

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <iomanip>

#include "Tool.h"

#include <TF1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TFitResult.h>
#include <TLegend.h>

using namespace std;

#define RESET   "\033[0m"
#define BLUE    "\033[34m"
#define GREEN   "\033[32m"

/**
 * \class GetDecayTime
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class GetDecayTime: public Tool {


    public:

        GetDecayTime(); ///< Simple constructor
        bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
        bool Execute(); ///< Executre function used to perform Tool perpose. 
        bool Finalise(); ///< Finalise funciton used to clean up resorces.

        vector<float> Decay_Ch0,Decay_Ch1,Decay_Ch2,Decay_Ch3,Decay_Ch4,Decay_Ch5,
                        Decay_Ch6,Decay_Ch7,Decay_Ch8,Decay_Ch9,Decay_Ch10,Decay_Ch11,
                        Decay_Ch12,Decay_Ch13,Decay_Ch14,Decay_Ch15;

        float guess_sigma, guess_t0, guess_n0, guess_tau0, guess_n1, guess_tau1, guess_n2, guess_tau2,guess_n3, guess_tau3, guess_offset,
                n_bins,bins_start,bins_end,fit_start,fit_end;    

        string fit_type;    

        TF1 *fit_function;
        TF1 *singleExpFunc0; // 4 parameters for Single_Exp
        TF1 *singleExpFunc1; // 4 parameters for Single_Exp
        TF1 *singleExpFunc2; // 4 parameters for Single_Exp
        TF1 *singleExpFunc3; // 4 parameters for Single_Exp


    private:
        int ROI_low;
        int ROI_high;

        float TriggerThresholdChannel;
        float ChannelThresholdChannel;
        float threshold_channel;
        float threshold_multiplier;
        int trg_ch;
        int words_in_data;
        int m_visibility;
        int fit_type_offset;
        int MultiplicityCut;
        bool log_style;

        void InitRoot();
        void FitDecay(int channel, vector<float> data);


};


#endif
