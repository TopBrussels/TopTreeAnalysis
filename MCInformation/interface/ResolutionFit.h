#ifndef ResolutionFit_h
#define ResolutionFit_h

#include <iostream>
#include <iomanip>
#include <vector>
#include "TH1F.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TLatex.h"

#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/Reconstruction/interface/FactorizedJetCorrector.h"
#include "TopTreeAnalysis/Reconstruction/interface/JetCorrectorParameters.h"

#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootMCParticle.h"

using namespace std;
using namespace TopTree;

// Calculation of the resolutions as done by Holger Enderle (slightly modified)
// see: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Bromo/TopAnalysis/TopUtils/bin/?pathrev=MAIN

class ResolutionFit {

  public:
    ResolutionFit();
    ResolutionFit(string label);
    ResolutionFit(const ResolutionFit &r);
    ~ResolutionFit();
    void Fill(TLorentzVector *lorentzVector, TRootMCParticle *mcParticle);
    void CalculateResolutions();
    void WritePlots(TFile *fout, bool savePNG = false, string pathPNG = string(""));
    void WriteResolutions(string file);
    void LoadResolutions(string file);
    
    float EtResolution(TLorentzVector* jet);
    float EtaResolution(TLorentzVector* jet);
    float ThetaResolution(TLorentzVector* jet);
    float PhiResolution(TLorentzVector* jet);
    float EtCorrection(TLorentzVector* jet);
    
  private:
    vector<double> ExtractSigmaMean(TH1* theHisto);

    FactorizedJetCorrector *jetCorr_;
    
    string label_;
    
    bool calculatedResolutions_;
    bool loadedResolutions_;
    
    Float_t towerBinning_[13]; // HCAL tower bins (for jets)
    Float_t jetPtBinning_[12];
    
    UInt_t nEtaBins_;
    UInt_t nPtBins_;
    
    TF1* bFuncSigma_[12];
    TF1* bFuncEtaSigma_[12];
    TF1* bFuncThetaSigma_[12];
    TF1* bFuncPhiSigma_[12];
    TF1* bFuncRelSigma_[12];
    
    TF1* bFuncMean_[12];
    TF1* bFuncRelMean_[12];
    
    TH1F* binCenterHisto_[143]; // size = (nEtaBins_*nPtBins_)+nPtBins_
    TH1F* resHisto_[143];
    TH1F* resHistoEta_[143];
    TH1F* resHistoTheta_[143];
    TH1F* resHistoPhi_[143];
    TH1F* resRelHisto_[143];
    TH1F* binCenterHistoIncl_[12];
    TH1F* resRelHistoIncl_[12];
        
    TCanvas* controlCan_[12];
    TCanvas* controlCanEta_[12];
    TCanvas* controlCanTheta_[12];
    TCanvas* controlCanPhi_[12];
    TCanvas* controlCanRel_[12];
    TCanvas* controlCanRelIncl_;

    TCanvas* bPtEtSigma_[12];
    TCanvas* bPtEtaSigma_[12];
    TCanvas* bPtThetaSigma_[12];
    TCanvas* bPtPhiSigma_[12];
    TCanvas* bPtEtRelSigma_[12];
    
    TCanvas* bPtEtMean_[12];
    TCanvas* bPtEtaMean_[12];
    TCanvas* bPtThetaMean_[12];
    TCanvas* bPtPhiMean_[12];
    TCanvas* bPtEtRelMean_[12];
    TCanvas* bPtEtRelMeanIncl_;
    
    TGraphErrors* bGraphSigma_[13];
    TGraphErrors* bGraphEtaSigma_[13];
    TGraphErrors* bGraphThetaSigma_[13];
    TGraphErrors* bGraphPhiSigma_[13];
    TGraphErrors* bGraphRelSigma_[13];
    
    TGraphErrors* bGraphMean_[13];
    TGraphErrors* bGraphEtaMean_[13];
    TGraphErrors* bGraphThetaMean_[13];
    TGraphErrors* bGraphPhiMean_[13];
    TGraphErrors* bGraphRelMean_[13];
    TGraphErrors* bGraphRelMeanIncl_;
};

#endif
