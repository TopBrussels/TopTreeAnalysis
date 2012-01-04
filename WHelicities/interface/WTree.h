#ifndef WTree_h
#define WTree_h

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "Rtypes.h"
#include "TObject.h"
#include "TVector3.h"
#include "TRef.h"
#include "TH2F.h"
#include "TLorentzVector.h"

using namespace std;

class WTree : public TObject
{
 public:
	
 WTree() :
  TObject()
    ,eventID_(0)
    ,runID_(0)
    ,lumiBlockID_(0)
    ,nPV_(0)
    ,nPUBXm1_(0)
    ,nPU_(0)
    ,nPUBXp1_(0)
    ,eventWeight_(0)
    ,chi2KinFit_()
    ,fittedLepton_()
    ,fittedNeutrino_()
    ,fittedLeptB_()
    ,fittedHadrB_()
    ,fittedLight1_()
    ,fittedLight2_()
    ,chi2FullKinFit_()
    ,fittedFullLepton_()
    ,fittedFullNeutrino_()
    ,fittedFullLeptB_()
    ,fittedFullHadrB_()
    ,fittedFullLight1_()
    ,fittedFullLight2_()
    ,chi2KinFitMassFit_()
    ,fittedLeptonMassFit_()
    ,fittedNeutrinoMassFit_()
    ,fittedLeptBMassFit_()
    ,fittedHadrBMassFit_()
    ,fittedLight1MassFit_()
    ,fittedLight2MassFit_()
    ,chi2FullKinFitMassFit_()
    ,fittedFullLeptonMassFit_()
    ,fittedFullNeutrinoMassFit_()
    ,fittedFullLeptBMassFit_()
    ,fittedFullHadrBMassFit_()
    ,fittedFullLight1MassFit_()
    ,fittedFullLight2MassFit_()
    ,hadrBJet_(0)
    ,hadrLJet1_(0)
    ,hadrLJet2_(0)
    ,leptBJet_(0)
    ,MET_()
    ,selectedJets_()
    ,bTagTCHE_()
    ,bTagTCHP_()
    ,bTagSSVHE_()
    ,bTagSSVHP_()
    ,bTagCSV_()
    ,muon_()
    ,hadrBQuark_()
    ,hadrLQuark1_()
    ,hadrLQuark2_()
    ,leptBQuark_()
    ,standardCosTheta_(-9999.)
    ,standardNeutrino_()
    ,standardLepton_()
    {;}
  
  ~WTree() {;}
  
  unsigned int eventID() const { return eventID_; }
  unsigned int runID() const { return runID_; }
  unsigned int lumiBlockID() const { return lumiBlockID_; }
  unsigned int nPV() const { return nPV_; }
  unsigned int nPUBXm1() const { return nPUBXm1_; }
  unsigned int nPU() const { return nPU_; }
  unsigned int nPUBXp1() const { return nPUBXp1_; }
  float eventWeight() const { return eventWeight_; }
  float chi2KinFit(int iCombi) { return chi2KinFit_[iCombi]; }
  TLorentzVector fittedLepton(int iCombi) {return fittedLepton_[iCombi];}
  TLorentzVector fittedNeutrino(int iCombi) {return fittedNeutrino_[iCombi];}
  TLorentzVector fittedLeptB(int iCombi) {return fittedLeptB_[iCombi];}
  TLorentzVector fittedHadrB(int iCombi) {return fittedHadrB_[iCombi];}
  TLorentzVector fittedLight1(int iCombi) {return fittedLight1_[iCombi];}
  TLorentzVector fittedLight2(int iCombi) {return fittedLight2_[iCombi];}
  float chi2FullKinFit(int iCombi) {return  chi2FullKinFit_[iCombi];}
  TLorentzVector fittedFullLepton(int iCombi) {return fittedFullLepton_[iCombi];}
  TLorentzVector fittedFullNeutrino(int iCombi) {return fittedFullNeutrino_[iCombi];}
  TLorentzVector fittedFullLeptB(int iCombi) {return fittedFullLeptB_[iCombi];}
  TLorentzVector fittedFullHadrB(int iCombi) {return fittedFullHadrB_[iCombi];}
  TLorentzVector fittedFullLight1(int iCombi) {return fittedFullLight1_[iCombi];}
  TLorentzVector fittedFullLight2(int iCombi) {return fittedFullLight2_[iCombi];}
  float chi2KinFitMassFit(int iCombi) {return chi2KinFitMassFit_[iCombi];}
  TLorentzVector fittedLeptonMassFit(int iCombi) {return fittedLeptonMassFit_[iCombi];}
  TLorentzVector fittedNeutrinoMassFit(int iCombi) {return fittedNeutrinoMassFit_[iCombi];}
  TLorentzVector fittedLeptBMassFit(int iCombi) {return fittedLeptBMassFit_[iCombi];}
  TLorentzVector fittedHadrBMassFit(int iCombi) {return fittedHadrBMassFit_[iCombi];}
  TLorentzVector fittedLight1MassFit(int iCombi) {return fittedLight1MassFit_[iCombi];}
  TLorentzVector fittedLight2MassFit(int iCombi) {return fittedLight2MassFit_[iCombi];}
  float chi2FullKinFitMassFit(int iCombi) {return chi2FullKinFitMassFit_[iCombi];}
  TLorentzVector fittedFullLeptonMassFit(int iCombi) {return fittedFullLeptonMassFit_[iCombi];}
  TLorentzVector fittedFullNeutrinoMassFit(int iCombi) {return fittedFullNeutrinoMassFit_[iCombi];}
  TLorentzVector fittedFullLeptBMassFit(int iCombi) {return fittedFullLeptBMassFit_[iCombi];}
  TLorentzVector fittedFullHadrBMassFit(int iCombi) {return fittedFullHadrBMassFit_[iCombi];}
  TLorentzVector fittedFullLight1MassFit(int iCombi) {return fittedFullLight1MassFit_[iCombi];}
  TLorentzVector fittedFullLight2MassFit(int iCombi) {return fittedFullLight2MassFit_[iCombi];}
  int hadrBJet() const { return hadrBJet_; } // index, according to MC
  int hadrLJet1() const { return hadrLJet1_; } // index, according to MC
  int hadrLJet2() const { return hadrLJet2_; } // index, according to MC
  int leptBJet() const { return leptBJet_; } // index, according to MC
  TLorentzVector met() const { return MET_; }
  vector<TLorentzVector> selectedJets() const { return selectedJets_; }
  TLorentzVector selectedJet(int i) const { return selectedJets_[i]; }
  vector<float> bTagTCHE() const { return bTagTCHE_; }
  vector<float> bTagTCHP() const { return bTagTCHP_; }
  vector<float> bTagSSVHE() const { return bTagSSVHE_; }
  vector<float> bTagSSVHP() const { return bTagSSVHP_; }
  vector<float> bTagCSV() const {return bTagCSV_;}
  TLorentzVector muon() const { return muon_; }
  TLorentzVector hadrBQuark() const { return hadrBQuark_; }
  TLorentzVector hadrLQuark1() const { return hadrLQuark1_; }
  TLorentzVector hadrLQuark2() const { return hadrLQuark2_; }
  TLorentzVector leptBQuark() const { return leptBQuark_; }
  float standardCosTheta() const {return standardCosTheta_;}
  TLorentzVector standardNeutrino() const {return standardNeutrino_;}
  TLorentzVector standardLepton() const {return standardLepton_;}
  
  void setEventID(unsigned int eventID) { eventID_ = eventID; }
  void setRunID(unsigned int runID) { runID_ = runID; }
  void setLumiBlockID(unsigned int lumiBlockID) { lumiBlockID_ = lumiBlockID; }
  void setNPV(unsigned int nPV) { nPV_ = nPV; }
  void setNPUBXm1(unsigned int nPUBXm1) { nPUBXm1_ = nPUBXm1; }
  void setNPU(unsigned int nPU) { nPU_ = nPU; }
  void setNPUBXp1(unsigned int nPUBXp1) { nPUBXp1_ = nPUBXp1; }
  void setEventWeight(float eventWeight) { eventWeight_ = eventWeight; }
  void setKinFitResults(vector<float> Chi2KinFitResults, vector<TLorentzVector> FittedLepton, vector<TLorentzVector> FittedNeutrino, vector<TLorentzVector> FittedLeptB, vector<TLorentzVector> FittedHadrB, vector<TLorentzVector> FittedLight1, vector<TLorentzVector> FittedLight2){
    for(unsigned int i=0; i<12; i++){
      chi2KinFit_[i] = Chi2KinFitResults[i];
      fittedLepton_[i] = FittedLepton[i];
      fittedNeutrino_[i] = FittedNeutrino[i];
      fittedLeptB_[i] = FittedLeptB[i];
      fittedHadrB_[i] = FittedHadrB[i];
      fittedLight1_[i] = FittedLight1[i];
      fittedLight2_[i] = FittedLight2[i];
    }
  }
  void setFullKinFitResults(vector<float> Chi2FullKinFitResults, vector<TLorentzVector> FittedFullLepton, vector<TLorentzVector> FittedFullNeutrino, vector<TLorentzVector> FittedFullLeptB, vector<TLorentzVector> FittedFullHadrB, vector<TLorentzVector> FittedFullLight1, vector<TLorentzVector> FittedFullLight2){
    for(unsigned int i=0; i<12; i++){
      chi2FullKinFit_[i] = Chi2FullKinFitResults[i];
      fittedFullLepton_[i] = FittedFullLepton[i];
      fittedFullNeutrino_[i] = FittedFullNeutrino[i];
      fittedFullLeptB_[i] = FittedFullLeptB[i];
      fittedFullHadrB_[i] = FittedFullHadrB[i];
      fittedFullLight1_[i] = FittedFullLight1[i];
      fittedFullLight2_[i] = FittedFullLight2[i];
    }
  }
 void setKinFitResultsMassFit(vector<float> Chi2KinFitResultsMassFit, vector<TLorentzVector> FittedLeptonMassFit, vector<TLorentzVector> FittedNeutrinoMassFit, vector<TLorentzVector> FittedLeptBMassFit, vector<TLorentzVector> FittedHadrBMassFit, vector<TLorentzVector> FittedLight1MassFit, vector<TLorentzVector> FittedLight2MassFit){
    for(unsigned int i=0; i<12; i++){
      chi2KinFitMassFit_[i] = Chi2KinFitResultsMassFit[i];
      fittedLeptonMassFit_[i] = FittedLeptonMassFit[i];
      fittedNeutrinoMassFit_[i] = FittedNeutrinoMassFit[i];
      fittedLeptBMassFit_[i] = FittedLeptBMassFit[i];
      fittedHadrBMassFit_[i] = FittedHadrBMassFit[i];
      fittedLight1MassFit_[i] = FittedLight1MassFit[i];
      fittedLight2MassFit_[i] = FittedLight2MassFit[i];
    }
  }
  void setFullKinFitResultsMassFit(vector<float> Chi2FullKinFitResultsMassFit, vector<TLorentzVector> FittedFullLeptonMassFit, vector<TLorentzVector> FittedFullNeutrinoMassFit, vector<TLorentzVector> FittedFullLeptBMassFit, vector<TLorentzVector> FittedFullHadrBMassFit, vector<TLorentzVector> FittedFullLight1MassFit, vector<TLorentzVector> FittedFullLight2MassFit){
    for(unsigned int i=0; i<12; i++){
      chi2FullKinFitMassFit_[i] = Chi2FullKinFitResultsMassFit[i];
      fittedFullLeptonMassFit_[i] = FittedFullLeptonMassFit[i];
      fittedFullNeutrinoMassFit_[i] = FittedFullNeutrinoMassFit[i];
      fittedFullLeptBMassFit_[i] = FittedFullLeptBMassFit[i];
      fittedFullHadrBMassFit_[i] = FittedFullHadrBMassFit[i];
      fittedFullLight1MassFit_[i] = FittedFullLight1MassFit[i];
      fittedFullLight2MassFit_[i] = FittedFullLight2MassFit[i];
    }
  }
  void setHadrBJet(int hadrBJet) { hadrBJet_ = hadrBJet; }
  void setHadrLJet1(int hadrLJet1) { hadrLJet1_ = hadrLJet1; }
  void setHadrLJet2(int hadrLJet2) { hadrLJet2_ = hadrLJet2; }
  void setLeptBJet(int leptBJet) { leptBJet_ = leptBJet; }
  void setMET(TLorentzVector MET) { MET_ = MET; }
  void setSelectedJets(vector<TLorentzVector> selectedJets) { selectedJets_ = selectedJets; }
  void setBTagTCHE(vector<float> bTagTCHE) { bTagTCHE_ = bTagTCHE; }
  void setBTagTCHP(vector<float> bTagTCHP) { bTagTCHP_ = bTagTCHP; }
  void setBTagSSVHE(vector<float> bTagSSVHE) { bTagSSVHE_ = bTagSSVHE; }
  void setBTagSSVHP(vector<float> bTagSSVHP) { bTagSSVHP_ = bTagSSVHP; }
  void setBTagCSV(vector<float> bTagCSV) {bTagCSV_ = bTagCSV;}
  void setMuon(TLorentzVector muon) { muon_ = muon; }
  void setHadrBQuark(TLorentzVector hadrBQuark) { hadrBQuark_ = hadrBQuark; }
  void setHadrLQuark1(TLorentzVector hadrLQuark1) { hadrLQuark1_ = hadrLQuark1; }
  void setHadrLQuark2(TLorentzVector hadrLQuark2) { hadrLQuark2_ = hadrLQuark2; }
  void setLeptBQuark(TLorentzVector leptBQuark) { leptBQuark_ = leptBQuark; }
  void setStandardCosTheta(float standardCosTheta) { standardCosTheta_ = standardCosTheta;}
  void setStandardNeutrino(TLorentzVector standardNeutrino) { standardNeutrino_ = standardNeutrino;}
  void setStandardLepton(TLorentzVector standardLepton) { standardLepton_ = standardLepton;}
  
 protected:
  
  unsigned int eventID_;
  unsigned int runID_;
  unsigned int lumiBlockID_;
  unsigned int nPV_;
  unsigned int nPUBXm1_;
  unsigned int nPU_;
  unsigned int nPUBXp1_;
  float eventWeight_;
  float chi2KinFit_[12];
  TLorentzVector fittedLepton_[12];
  TLorentzVector fittedNeutrino_[12];
  TLorentzVector fittedLeptB_[12];
  TLorentzVector fittedHadrB_[12];
  TLorentzVector fittedLight1_[12];
  TLorentzVector fittedLight2_[12];
  float chi2FullKinFit_[12];
  TLorentzVector fittedFullLepton_[12];
  TLorentzVector fittedFullNeutrino_[12];
  TLorentzVector fittedFullLeptB_[12];
  TLorentzVector fittedFullHadrB_[12];
  TLorentzVector fittedFullLight1_[12];
  TLorentzVector fittedFullLight2_[12];
  float chi2KinFitMassFit_[12];
  TLorentzVector fittedLeptonMassFit_[12];
  TLorentzVector fittedNeutrinoMassFit_[12];
  TLorentzVector fittedLeptBMassFit_[12];
  TLorentzVector fittedHadrBMassFit_[12];
  TLorentzVector fittedLight1MassFit_[12];
  TLorentzVector fittedLight2MassFit_[12];  
  float chi2FullKinFitMassFit_[12];
  TLorentzVector fittedFullLeptonMassFit_[12];
  TLorentzVector fittedFullNeutrinoMassFit_[12];
  TLorentzVector fittedFullLeptBMassFit_[12];
  TLorentzVector fittedFullHadrBMassFit_[12];
  TLorentzVector fittedFullLight1MassFit_[12];
  TLorentzVector fittedFullLight2MassFit_[12];  
  int hadrBJet_; //index according to MC
  int hadrLJet1_;
  int hadrLJet2_;
  int leptBJet_;
  TLorentzVector MET_;
  vector<TLorentzVector> selectedJets_; // all selected jets
  vector<float> bTagTCHE_; // indices like selectedJets indices
  vector<float> bTagTCHP_;
  vector<float> bTagSSVHE_;
  vector<float> bTagSSVHP_;
  vector<float> bTagCSV_;
  TLorentzVector muon_;
  TLorentzVector hadrBQuark_;
  TLorentzVector hadrLQuark1_;
  TLorentzVector hadrLQuark2_;
  TLorentzVector leptBQuark_;
  float standardCosTheta_;
  TLorentzVector standardNeutrino_;
  TLorentzVector standardLepton_;
  
  ClassDef (WTree,2);
};

#endif

