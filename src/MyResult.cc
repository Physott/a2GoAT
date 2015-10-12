#include "MyResult.h"



MyResult::MyResult(const char *_Name, TFile *input, TFile *output, const int _Color)  :
    in(input),
    out(output),
    result(_Name, out, _Color)
{
}

MyResult::~MyResult()
{

}

Bool_t	MyResult::Start(const int summedTaggerChannels, const int summedIMBins)
{
    if(!in)     return false;
    if(!in->IsOpen())   return false;
    if(!out)    return false;
    if(!out->IsOpen())   return false;
    if(!out->IsWritable())   return false;

    result.SetFile(in, summedTaggerChannels);
    result.RebinIM(summedIMBins);
    result.Draw();

	return kTRUE;
}

