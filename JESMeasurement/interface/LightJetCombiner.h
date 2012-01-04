#ifndef LightJetCombiner_h
#define LightJetCombiner_h

// Class that takes care of the Jet-parton matching and of the MVA-based jet combination

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdio>

#include "TLorentzVector.h"
#include "TH1F.h"

#include <sys/stat.h>

#include "TopTreeAnalysis/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysis/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysis/JESMeasurement/interface/Monster.h"
#include "TopTreeAnalysis/JESMeasurement/interface/MonsterTools.h"
#include "TopTreeAnalysis/JESMeasurement/interface/FullKinFit.h"

using namespace std;

struct MVAValues {
  string MVAAlgo;
  float MVAValue;
  int WJet1;
  int WJet2;
  int HadrBJet;
  int LeptBJet;
} ;

class LightJetCombiner {

  public:
    LightJetCombiner(bool trainMVA, float Luminosity, const vector<Dataset*>& datasets, bool measureTopMass = true, string MVAMethod = "Likelihood");
    ~LightJetCombiner();
    void ProcessEvent(Dataset* dataSet, Monster* monster, float scaleFactor=1);

    pair<float, vector<unsigned int> > getMVAValue(string MVAMethod, int rank); // rank 1 means the highest MVA value for this method, rank 2 the second highest, ...
    void Write(TFile* fout, bool savePNG = false, string pathPNG = string(""));
    
  private:
    map<string, MultiSamplePlot*> MSPlot_;
    map<string, TH1F*> histo1D_;
    map<string, TH2F*> histo2D_;
    map<string,TGraphAsymmErrors*> graphAsymmErr_;
    
    vector<Dataset*> datasets_;
    float Luminosity_;
    
    bool trainMVA_; // true if the MVA needs to be trained, false if the trained MVA will be used to compute stuff
    MVATrainer* trainer_;
    MVAComputer* computer_;

    vector<MVAValues> vectorMVA_;
    
    bool measureTopMass_;
    bool EndJobWasRun_;
    
    TH2F* dummyMonster_;
};

#endif
