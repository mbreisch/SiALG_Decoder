#ifndef LoadBinaryFile_H
#define LoadBinaryFile_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "Tool.h"

using namespace std;

#define MAX_NUM_CHANNELS 16

/**
 * \class LoadBinaryFile
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/

class LoadBinaryFile: public Tool 
{
    public:

        LoadBinaryFile(); ///< Simple constructor
        bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
        bool Execute(); ///< Executre function used to perform Tool perpose. 
        bool Finalise(); ///< Finalise funciton used to clean up resorces.


        std::string File_Prefix;
        std::string File_Sufix;
        std::string FullFilePath;
        std::string Trigger_0;
        std::string Trigger_1;

    private:

        int words_in_header;
        int bytes_in_headerword;
        int words_in_data;
        int bytes_in_dataword;

        vector<unsigned int> tmp_H;
        vector<float> tmp_D;

        vector<unsigned int> LoadHeader(int i_channel);
        vector<float> LoadData(int i_channel);

};
#endif
