#include "Factory.h"

Tool* Factory(std::string tool) {
Tool* ret=0;

// if (tool=="Type") tool=new Type;
if (tool=="DummyTool") ret=new DummyTool;
if (tool=="LoadBinaryFile") ret=new LoadBinaryFile;
if (tool=="StoreEvent") ret=new StoreEvent;
if (tool=="BaselineSubstraction") ret=new BaselineSubstraction;
  if (tool=="PulseIdentifier") ret=new PulseIdentifier;
return ret;
}
