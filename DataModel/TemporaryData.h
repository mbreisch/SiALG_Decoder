#ifndef TEMPORARYDATA_H
#define TEMPORARYDATA_H

#include <SerialisableObject.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>

#include "TTree.h"
#include "TFile.h"

using namespace std;

class TemporaryData{

    friend class boost::serialization::access;

 public:

    TemporaryData();
    ~TemporaryData();

    std::string Path_In;
    std::string Path_Out;

    map<int,std::ifstream> FileMap;
    map<int,vector<unsigned int>> ParsedMap_Header;
    map<int,vector<float>> ParsedMap_Data;
    map<int,vector<int>> PeakPositions;
    map<int,vector<float>> PeakMinima;
    map<int,vector<float>> TTS_Map;

    vector<int> ListOfChannels;

    int EventCounter;
    int Runtime;
    bool Stop;

    TFile* RootFile_Event;
    TTree *TTree_Event;

    TFile* RootFile_Analysis;
    TTree *TTree_Analysis_PHD;
    TTree *TTree_Analysis_TTS;

 private:
 
};

#endif
