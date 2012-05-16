#ifndef TopFCNC_KinFit_h
#define TopFCNC_KinFit_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <sys/stat.h>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TH2F.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TopTreeAnalysis/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysis/KinFitter/interface/TKinFitter.h"
#include "TopTreeAnalysis/KinFitter/interface/TFitConstraintM.h"
#include "TopTreeAnalysis/KinFitter/interface/TFitParticleEtThetaPhiEMomFix.h"
#include "TopTreeAnalysis/Reconstruction/interface/FactorizedJetCorrector.h"
#include "TopTreeAnalysis/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/TopFCNC/interface/TopFCNC_Evt.h"

#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootMuon.h"

using namespace std;
using namespace TopTree;

class TopFCNC_KinFit {

  public:
    TopFCNC_KinFit(Dataset* dataset, ResolutionFit *resFitLeptons, ResolutionFit *resFitLightJets, ResolutionFit *resFitBJets);
    ~TopFCNC_KinFit();
    void FitEvent(TopFCNC_Evt *topFCNC_Evt, float WMass = 80.4, float Zmass = 91.2, float topMass = 172.5);
    void Write(TFile* fout, bool savePNG = false, string pathPNG = string(""));

    Double_t GetProb() {return Prob_;};
    Double_t GetChi2() {return Chi2_;};
    Double_t GetNdof() {return Ndof_;};
    
    void SetMaxNbIter(Int_t MaxNbIter)    {MaxNbIter_ = MaxNbIter;}
    void SetMaxDeltaS(Double_t MaxDeltaS) {MaxDeltaS_ = MaxDeltaS;}
    void SetMaxF(Double_t MaxF)           {MaxF_ = MaxF;}

    void SetVerbosity(Bool_t verbose) {verbose_ = verbose;};
    void SetFitVerbosity(Bool_t fit_verbose) {fit_verbose_ = fit_verbose;};
  
  private:
    map<string, TH1F*> histo1D_;
    Dataset       *dataset_;

    ResolutionFit *resFitLeptons_;
    ResolutionFit *resFitLightJets_;
    ResolutionFit *resFitBJets_;

    Double_t       Prob_;
    Double_t       Chi2_;
    Int_t          Ndof_;
    Int_t          MaxNbIter_;  // Maximum number of iterations
    Double_t       MaxDeltaS_;  // Convergence criterium for deltaS
    Double_t       MaxF_;       // Convergence criterium for F
    
    Bool_t         verbose_;
    Bool_t         fit_verbose_;
};

#endif
