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

    vector<float> PusleHeight_Ch0,PusleHeight_Ch1,PusleHeight_Ch2,PusleHeight_Ch3,PusleHeight_Ch4,PusleHeight_Ch5,
                    PusleHeight_Ch6,PusleHeight_Ch7,PusleHeight_Ch8,PusleHeight_Ch9,PusleHeight_Ch10,PusleHeight_Ch11,
                    PusleHeight_Ch12,PusleHeight_Ch13,PusleHeight_Ch14,PusleHeight_Ch15;

    vector<float> GetSlicedDataFromROI(vector<float> data, int startIndex, int endIndex);

    private:
        int ROI_low;
        int ROI_high;


};
#endif
