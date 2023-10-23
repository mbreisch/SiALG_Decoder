#include "Factory.h"

Tool* Factory(std::string tool) {
Tool* ret=0;

// if (tool=="Type") tool=new Type;
if (tool=="BaselineSubstraction") ret=new BaselineSubstraction;
if (tool=="DummyTool") ret=new DummyTool;
if (tool=="EndChain") ret=new EndChain;
if (tool=="GetCharge") ret=new GetCharge;
if (tool=="GetDecayTime") ret=new GetDecayTime;
if (tool=="GetPhDistribution") ret=new GetPhDistribution;
if (tool=="GetTTS") ret=new GetTTS;
if (tool=="LoadBinaryFile") ret=new LoadBinaryFile;
if (tool=="PulseIdentifier") ret=new PulseIdentifier;
if (tool=="StoreEvent") ret=new StoreEvent;
  if (tool=="LoadRootFile") ret=new LoadRootFile;
return ret;
}
