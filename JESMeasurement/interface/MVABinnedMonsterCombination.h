#ifndef MVABinnedMonsterCombination_h
#define MVABinnedMonsterCombination_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include "TROOT.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TFile.h"

#include "TopTreeAnalysis/JESMeasurement/interface/MonsterCombination.h"
#include "TopTreeAnalysis/JESMeasurement/interface/JESTools.h"
#include "TopTreeAnalysis/JESMeasurement/interface/Monster.h"

using namespace std;

class MVABinnedMonsterCombination
{
  public:
    MVABinnedMonsterCombination(string title, bool topMass, vector<float> binning);
    ~MVABinnedMonsterCombination();
    void Fill(Monster* monster, float weight=1, int jetCombi=0);
    void FillExpected(Monster* monster, float binValue, float weight=1);
    void Write(TFile* f, string pathforPNG, bool templateFit, TH1F* goodMVAtemplate, TH1F* badMVAtemplate, string xTitle = "G / ( G + B )", bool silent = true);
  
  private:
    string title_;
    bool topMass_;
    
    vector<float> binning_;
    vector<MonsterCombination*> monsterCombis_;
    vector<TH1F*> xBinHistos_;
    
    map<string, vector<TH1F*> > mvaBinnedVars_;
    
    TH1F* maxMVA_;
    TH1F* expDElHisto_;
    TH1F* expDEbHisto_;
};

#endif
