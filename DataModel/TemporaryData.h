#ifndef TEMPORARYDATA_H
#define TEMPORARYDATA_H

#include <SerialisableObject.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>

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
    map<int,vector<float>> Charge_Map;
    map<int,vector<float>> Decay_Map;

    vector<int> ListOfChannels;

    int EventCounter;
    int Runtime;
    bool Stop;
    bool EndOfRun;
    bool NewRun;

    TFile* RootFile_Event;
    TTree *TTree_Event;

    TFile* RootFile_Analysis;
    TTree *TTree_Analysis_PHD;
    TTree *TTree_Analysis_TTS;
    TTree *TTree_Analysis_Charge;
    TTree *TTree_Analysis_DecayTime;
    
 private:
 
};

#endif
