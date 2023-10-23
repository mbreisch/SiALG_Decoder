#ifndef LoadRootFile_H
#define LoadRootFile_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include <TFile.h>
#include <TTree.h>

#include "Tool.h"

using namespace std;

#define MAX_NUM_CHANNEL 16


/**
 * \class LoadRootFile
 *
 * This is a balnk template for a Tool used by the script to generate a new custom tool. Please fill out the descripton and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
*/

class LoadRootFile: public Tool 
{
    public:

        LoadRootFile(); ///< Simple constructor
        bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
        bool Execute(); ///< Executre function used to perform Tool perpose. 
        bool Finalise(); ///< Finalise funciton used to clean up resorces.


    private:

        std::ifstream file;
        std::string Treename;
        std::string Branchname;
        std::string line;

};


#endif
