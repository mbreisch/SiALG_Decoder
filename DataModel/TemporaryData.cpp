#include <TemporaryData.h>

TemporaryData::TemporaryData()
{
    Stop= false;
    EndOfRun = false;
    NewRun = true;
    ListOfChannels = {};
}

TemporaryData::~TemporaryData()
{}