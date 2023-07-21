#ifndef GetCharge_H
#define GetCharge_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TTree.h"
#include "TFile.h"
#include <TMath.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TApplication.h>

using namespace std;
/**
 * \class GetCharge
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/

class GetCharge: public Tool 
{
    public:

    GetCharge(); ///< Simple constructor
    bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
    bool Execute(); ///< Executre function used to perform Tool perpose. 
    bool Finalise(); ///< Finalise funciton used to clean up resorces.

    vector<float> Charge_Ch0,Charge_Ch1,Charge_Ch2,Charge_Ch3,Charge_Ch4,Charge_Ch5,
                    Charge_Ch6,Charge_Ch7,Charge_Ch8,Charge_Ch9,Charge_Ch10,Charge_Ch11,
                    Charge_Ch12,Charge_Ch13,Charge_Ch14,Charge_Ch15;

    vector<float> GetSlicedDataFromROI(vector<float> data, int startIndex, int endIndex);

    private:
        int ROI_low;
        int ROI_high;
        int vis;
        int ECMAX;

};
#endif
