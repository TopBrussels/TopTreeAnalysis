#ifndef MonsterCombination_h
#define MonsterCombination_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include "TROOT.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"

#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/JESMeasurement/interface/JESTools.h"
#include "TopTreeAnalysis/JESMeasurement/interface/MonsterTools.h"
#include "TopTreeAnalysis/JESMeasurement/interface/Monster.h"
#include "TopTreeAnalysis/JESMeasurement/interface/FullKinFit.h"

using namespace std;

class MonsterCombination
{
  public:
    MonsterCombination(string title, bool topMass);
    ~MonsterCombination();
    void addMonster(Monster* monster, float weight=1, int jetCombi=0);
    void addExpected(Monster* monster, float weight=1, int jetCombi=0);
    TH2F* AverageMonster();
    void SubTractAvMonster(TH2F* avMonster, float nBadCombis);
    pair<float,std::pair<float,float> > GetExpectedDEl(); //get the expected DEl corrections
    pair<float,std::pair<float,float> > GetExpectedDEb(); //get the expected DEb corrections
    vector< vector<float> > GetEstimation(string pathforPNG, TFile* f, bool silent = true, string calibCurveName = ""); //get the estimated corrections
    
  private:
    TH2F* superMonster_;
    string title_;
    bool topMass_; // measuring the top mass or the JES
    float nEventsWeighted_;
    
    TH1F expDElHisto_;
    TH1F expDEbHisto_;
    
    vector<TLorentzVector> hadrLJet1_;
    vector<TLorentzVector> hadrLJet2_;
    vector<TLorentzVector> hadrBJet_;
};

#endif
