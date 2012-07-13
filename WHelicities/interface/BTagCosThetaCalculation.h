#ifndef BTagCosThetaCalculation_h
#define BTagCosThetaCalculation_h

#include "TLorentzVector.h"
#include <Math/VectorUtil.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../../Selection/interface/SelectionTable.h"
#include "../../Tools/interface/PlottingTools.h"
#include "../../Tools/interface/MultiSamplePlot.h"
#include "../../Tools/interface/TTreeLoader.h"
#include "../../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../../Content/interface/AnalysisEnvironment.h"
#include "../../Content/interface/Dataset.h"
#include "../../MCInformation/interface/JetPartonMatching.h"
#include "../../Tools/interface/JetTools.h"
#include "TVector3.h"
#include "TH1.h"
#include "Riostream.h"
#include "../../Reconstruction/interface/JetCorrectorParameters.h"
#include "../../Reconstruction/interface/JetCorrectionUncertainty.h"

using namespace std;

class BTagCosThetaCalculation{
  
 public:
  BTagCosThetaCalculation();
  ~BTagCosThetaCalculation();
  
  float Calculation(TLorentzVector muon, TLorentzVector Neutrino, TLorentzVector leptonicBJet);
  float CalcOrigKins(int BLeptonicIndex, int BHadronicIndex,  TLorentzVector muon, vector<TLorentzVector> selectedJets, float MassW, float MassTop);
  
  TLorentzVector GetNeutrino() {return Neutrino;}
  
 private:
  TLorentzVector Neutrino;
  
};

#endif

