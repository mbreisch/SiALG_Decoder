#ifndef TimeStampTest_H
#define TimeStampTest_H

#include <string>
#include <iostream>
#include <map>
#include <vector>

#include "Tool.h"

#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>

/**
 * \class TimeStampTest
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class TimeStampTest: public Tool 
{
    public:

        TimeStampTest(); ///< Simple constructor
        bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
        bool Execute(); ///< Executre function used to perform Tool perpose. 
        bool Finalise(); ///< Finalise funciton used to clean up resorces.

        vector<unsigned int> Timestamp_Ch0,Timestamp_Ch1,Timestamp_Ch2,Timestamp_Ch3,Timestamp_Ch4,Timestamp_Ch5,
                        Timestamp_Ch6,Timestamp_Ch7,Timestamp_Ch8,Timestamp_Ch9,Timestamp_Ch10,Timestamp_Ch11,
                        Timestamp_Ch12,Timestamp_Ch13,Timestamp_Ch14,Timestamp_Ch15;

    private:
        unsigned int TriggerThresholdChannel;
        unsigned int ChannelThresholdChannel;
        int trg_ch;

        void InitRoot();

};

#endif
