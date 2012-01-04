#ifndef WReconstruction_h
#define WReconstruction_h

#include "TMatrixD.h"

//Use KinFitter present in CMSSW:
//#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
//#include "PhysicsTools/KinFitter/interface/TFitParticleEtThetaPhi.h"
//#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"

#include "TopTreeAnalysis/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysis/KinFitter/interface/TKinFitter.h"
#include "TopTreeAnalysis/KinFitter/interface/TFitParticleEtThetaPhi.h"
#include "TopTreeAnalysis/KinFitter/interface/TFitConstraintM.h"

// system include files
#include <memory>
#include <vector>
#include <iostream>

using namespace std;

class WReconstruction{

 public:
  WReconstruction();
  ~WReconstruction();

  void Initializing(bool Simulation,std::vector<float> JetsPx, std::vector<float> JetsPy, std::vector<float> JetsPz, std::vector<float> JetsE,std::vector<float> TCHEbTagValues,std::vector<float> SSVHEbTagValues,std::vector<float> EtResolutionLightJets,std::vector<float> ThetaResolutionLightJets,std::vector<float> PhiResolutionLightJets,std::vector<float> EtResolutionBJets,std::vector<float> ThetaResolutionBJets,std::vector<float> PhiResolutionBJets);
  std::vector<int> BTagAnalysis();
  float KinFitCalculation(int jet1Number, int jet2Number, int bhadrNumber);
  int CalculateLowestChiSquared(int sizeLoop);
  std::vector<int> QuarkDistribution(int CombinationValue);
  bool ControlChecks();
  float applicationBTag(int bJet1, int bJet2);
  
 private:
  std::vector<int> CorrectJets;

  std::vector<float> jetsPx_;
  std::vector<float> jetsPy_;
  std::vector<float> jetsPz_;
  std::vector<float> jetsE_;
  std::vector<float> muonsPx_;
  std::vector<float> muonsPy_;
  std::vector<float> muonsPz_;
  std::vector<float> muonsE_;
  std::vector<float> TCHEbTags_;
  std::vector<float> SSVHEbTags_;
  std::vector<float> EtResLight_;
  std::vector<float> ThetaResLight_;
  std::vector<float> PhiResLight_;
  std::vector<float> EtResB_;
  std::vector<float> ThetaResB_;
  std::vector<float> PhiResB_;

  std::vector<int> JetNumbers;
  std::vector<int> JetsEvent;

  int JetSize;
  int TCHEbTagSize;
  int SSVHEbTagSize;

  //  ResolutionFit *resFitLightJets_;
  //  ResolutionFit *resFitBJets_;  
  float TopMassKinFit;
  float WMassKinFit;

  float ChiSquared[12];
  float ChiSquaredValue;
  int NumberCombinations;
  std::vector<float> bTagBoolean;
  float bTagSelected;
  double btag[12];
  bool ProblemEvent;
  int QuarkOneIndex[12];
  int QuarkTwoIndex[12];
  int BHadronicIndex[12];
  int BLeptonicIndex;      
  double BTagValueMVA;
  int BLeptonicIndexMVA[12];
  int MVACombination[2];
  int UsedCombination;
  int MVAComb;
  float ChiSquaredFit;
  int CombinationOfMVA;

  float MassW;
  float MassTop;
  float SigmaW;  
  float SigmaTop;
};

#endif
