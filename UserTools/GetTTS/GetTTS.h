#ifndef GetTTS_H
#define GetTTS_H

#include <string>
#include <iostream>

#include "Tool.h"

#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>


/**
 * \class GetTTS
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/

class GetTTS: public Tool 
{
    public:

        GetTTS(); ///< Simple constructor
        bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
        bool Execute(); ///< Executre function used to perform Tool perpose. 
        bool Finalise(); ///< Finalise funciton used to clean up resorces.

        vector<float> TTS_Ch0,TTS_Ch1,TTS_Ch2,TTS_Ch3,TTS_Ch4,TTS_Ch5,
                        TTS_Ch6,TTS_Ch7,TTS_Ch8,TTS_Ch9,TTS_Ch10,TTS_Ch11,
                        TTS_Ch12,TTS_Ch13,TTS_Ch14,TTS_Ch15;

        vector<float> GetSlice(vector<float> data, int startIndex, int endIndex);
        float GetTriggerThreshold(int trg_ch);
        float GetChannelThreshold(vector<float> data, int position,float amplitude);

    private:
        int ROI_low;
        int ROI_high;

        float TriggerThresholdChannel;
        float ChannelThresholdChannel;
        float threshold_channel;
        float threshold_multiplier;
        int trg_ch;
        int words_in_data;
        int MultiplicityCut;

        void InitRoot();

};
#endif
