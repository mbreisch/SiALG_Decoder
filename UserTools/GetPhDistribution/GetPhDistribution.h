#ifndef GetPhDistribution_H
#define GetPhDistribution_H

#include <string>
#include <iostream>
#include <map>
#include <vector>

#include "Tool.h"

#include "TTree.h"
#include "TFile.h"
#include <TMath.h>

using namespace std;

#define MAX_NUM_CHANNEL 16

/**
 * \class GetPhDistribution
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/

class GetPhDistribution: public Tool 
{
    public:

    GetPhDistribution(); ///< Simple constructor
    bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
    bool Execute(); ///< Executre function used to perform Tool perpose. 
    bool Finalise(); ///< Finalise funciton used to clean up resorces.

    vector<float> PulseHeight_Ch0,PulseHeight_Ch1,PulseHeight_Ch2,PulseHeight_Ch3,PulseHeight_Ch4,PulseHeight_Ch5,
                    PulseHeight_Ch6,PulseHeight_Ch7,PulseHeight_Ch8,PulseHeight_Ch9,PulseHeight_Ch10,PulseHeight_Ch11,
                    PulseHeight_Ch12,PulseHeight_Ch13,PulseHeight_Ch14,PulseHeight_Ch15;

    vector<float> GetSlicedDataFromROI(vector<float> data, int startIndex, int endIndex);

    private:
        int ROI_low;
        int ROI_high;

        void InitRoot();

};
#endif
