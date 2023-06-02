#ifndef StoreEvent_H
#define StoreEvent_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TTree.h"
#include "TFile.h"


using namespace std;

/**
 * \class StoreEvent
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/

class StoreEvent: public Tool 
{
    public:

    StoreEvent(); ///< Simple constructor
    bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
    bool Execute(); ///< Executre function used to perform Tool perpose. 
    bool Finalise(); ///< Finalise funciton used to clean up resorces.

    std::string Path_Out;

    std::vector<unsigned int> Header_ch0,Header_ch1,Header_ch2,Header_ch3,Header_ch4,Header_ch5,Header_ch6,Header_ch7,Header_ch8,
        Header_ch9,Header_ch10,Header_ch11,Header_ch12,Header_ch13,Header_ch14,Header_ch15; 
    std::vector<float> Data_ch0,Data_ch1,Data_ch2,Data_ch3,Data_ch4,Data_ch5,Data_ch6,Data_ch7,Data_ch8,
        Data_ch9,Data_ch10,Data_ch11,Data_ch12,Data_ch13,Data_ch14,Data_ch15; 
    std::vector<unsigned int> Trigger_header_g0,Trigger_header_g1;
    std::vector<float> Trigger_data_g0,Trigger_data_g1;
    int EventID;

    void printMap(const std::map<int, std::vector<std::vector<unsigned int>>>& myMap);

    private:

};
#endif
