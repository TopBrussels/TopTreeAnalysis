#ifndef ExpCorrCalculator_h
#define ExpCorrCalculator_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <sys/stat.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TVector3.h"

#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/JESMeasurement/interface/JESTools.h"

#include "TopTreeProducer/interface/TRootJet.h"

using namespace std;

class ExpCorrCalculator {

  public:
    ExpCorrCalculator(string name);
    ~ExpCorrCalculator();
    void FillLightJet(float DEl, float PtLight = 0, float PtB = 0);
    void FillBJet(float DEb, float PtLight = 0, float PtB = 0);
    void Fill3Jets(vector<TLorentzVector> TLV);
    void Fill3MCParticles(vector<TLorentzVector> TLV);
    void FillMuons(vector<TLorentzVector> TLV);
    void FillElectrons(vector<TLorentzVector> TLV);
    void FillMaxMVA(float maxMVA);
    pair<float,float> GetFittedMass(TH1F* histo);
    void Calculate(bool printToScreen = false);
    void Write(TFile* fout, bool savePNG = false, string pathPNG = string(""));
    pair<float,float> GetExpCorrL();
    pair<float,float> GetExpCorrB();
  
  private:
    map<string, TH1F*> histo1D_;
    map<string, TH2F*> histo2D_;
    map<string,TProfile*> profile_;
    map<string,TGraphErrors*> graphErr_;
    
    string name_;
    
    Float_t jetPtBinning_[5];
    UInt_t nPtBins_;
    TF1* diffFunction_;
    TF1* dEbFunction_;
    
    bool wasCalculated_;
    
    pair<float,float> inclLCorr_;
    pair<float,float> inclBCorr_;
    
    vector<float> maxMVAs_;
    
    vector< vector<TLorentzVector> > jetTLVs_;
    vector< vector<TLorentzVector> > mcParticlesTLVs_;
    vector< vector<TLorentzVector> > muonTLVs_;
    vector< vector<TLorentzVector> > electronTLVs_;
};

#endif
