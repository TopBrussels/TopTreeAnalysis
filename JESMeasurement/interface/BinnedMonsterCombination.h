#ifndef BinnedMonsterCombination_h
#define BinnedMonsterCombination_h

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

#include "TopTreeAnalysis/JESMeasurement/interface/MVABinnedMonsterCombination.h"
#include "TopTreeAnalysis/JESMeasurement/interface/MonsterCombination.h"
#include "TopTreeAnalysis/JESMeasurement/interface/JESTools.h"
#include "TopTreeAnalysis/JESMeasurement/interface/Monster.h"

using namespace std;

class BinnedMonsterCombination
{
  public:
    BinnedMonsterCombination(string title, bool topMass, vector<float> binning, bool doMvaBinned=false);
    ~BinnedMonsterCombination();
    void Fill(Monster* monster, float binValue, float weight=1, int jetCombi=0);
    void FillExpected(Monster* monster, float binValue, float weight=1);
    TH2F* AverageMonster(float binValue);
    void SubTractAvMonster(TH2F* avMonster, float binValue, float nBadCombis);
    void Write(TFile* f, string pathforPNG, string xTitle = "", bool silent = true);
  
  private:
    string title_;
    bool topMass_;
    bool doMvaBinned_;
    
    map<string, vector<TH1F*> > binnedVars_;
    
    vector<float> binning_;
    vector<MonsterCombination*> monsterCombis_;
    vector<MVABinnedMonsterCombination*> monsterCombisMVAbinned_;
    vector<TH1F*> xBinHistos_;
};

#endif

