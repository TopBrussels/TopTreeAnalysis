#ifndef BTagCosThetaCalculation_h
#define BTagCosThetaCalculation_h

#include "TLorentzVector.h"
#include <Math/VectorUtil.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TVector3.h"
#include "TH1.h"
#include "Riostream.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"

#include <string>
#include <vector>
#include <iostream>

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

