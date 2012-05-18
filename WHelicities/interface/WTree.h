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
    //    ,eventWeight_(0)
    ,hadrKinFitChiSq_()
    ,hadrKinFitHadrB_()
    ,hadrKinFitLight1_()
    ,hadrKinFitLight2_()
    ,hadrAndLeptWOnlyKinFitChiSq_()
    ,hadrAndLeptWOnlyKinFitLepton_()
    ,hadrAndLeptWOnlyKinFitNeutrino_()
    ,hadrAndLeptWOnlyKinFitHadrB_()
    ,hadrAndLeptWOnlyKinFitLight1_()
    ,hadrAndLeptWOnlyKinFitLight2_()
    ,hadrAndLeptKinFitChiSq_()
    ,hadrAndLeptKinFitLepton_()
    ,hadrAndLeptKinFitNeutrino_()
    ,hadrAndLeptKinFitLeptB_()
    ,hadrAndLeptKinFitHadrB_()
    ,hadrAndLeptKinFitLight1_()
    ,hadrAndLeptKinFitLight2_()
    ,hadrBJet_(0)
    ,hadrLJet1_(0)
    ,hadrLJet2_(0)
    ,leptBJet_(0)
    ,MET_()
    ,selectedJets_()
    ,lepton_()
    ,bTagTCHE_()
    ,bTagTCHP_()
    ,bTagSSVHE_()
    ,bTagSSVHP_()
    ,bTagCSV_()
    ,bTagJP_()
    ,bTagJBP_()
    ,hadrBQuark_()
    ,hadrLQuark1_()
    ,hadrLQuark2_()
    ,leptBQuark_()
    ,standardCosTheta_(-9999.)
    ,standardNeutrino_()
    ,standardLepton_()
    ,recMassW_()
    ,recMassTop_()
    ,recMassWLept_()
    ,recMassTopLept_()
    ,relIso_()
    ,isoMuTriggerBool_()
    ,muonTriggerValue_()
    ,electronTriggerValue_()
    ,deltaRJetLepton_()
    ,deltaRMuonJet_()
    ,deltaRElectronJet_()
    {;}
  
  ~WTree() {;}
  
  unsigned int eventID() const { return eventID_; }
  unsigned int runID() const { return runID_; }
  unsigned int lumiBlockID() const { return lumiBlockID_; }
  unsigned int nPV() const { return nPV_; }
  unsigned int nPUBXm1() const { return nPUBXm1_; }
  unsigned int nPU() const { return nPU_; }
  unsigned int nPUBXp1() const { return nPUBXp1_; }
  //  float eventWeight() const { return eventWeight_; }
  float hadrKinFitChiSq(int iCombi) { return hadrKinFitChiSq_[iCombi];}
  TLorentzVector hadrKinFitHadrB(int iCombi) {return hadrKinFitHadrB_[iCombi];}
  TLorentzVector hadrKinFitLight1(int iCombi) {return hadrKinFitLight1_[iCombi];}
  TLorentzVector hadrKinFitLight2(int iCombi) {return hadrKinFitLight2_[iCombi];}
  float hadrAndLeptWOnlyKinFitChiSq(int iCombi) { return hadrAndLeptWOnlyKinFitChiSq_[iCombi]; }
  TLorentzVector hadrAndLeptWOnlyKinFitLepton(int iCombi) {return hadrAndLeptWOnlyKinFitLepton_[iCombi];}
  TLorentzVector hadrAndLeptWOnlyKinFitNeutrino(int iCombi) {return hadrAndLeptWOnlyKinFitNeutrino_[iCombi];}
  TLorentzVector hadrAndLeptWOnlyKinFitHadrB(int iCombi) {return hadrAndLeptWOnlyKinFitHadrB_[iCombi];}
  TLorentzVector hadrAndLeptWOnlyKinFitLight1(int iCombi) {return hadrAndLeptWOnlyKinFitLight1_[iCombi];}
  TLorentzVector hadrAndLeptWOnlyKinFitLight2(int iCombi) {return hadrAndLeptWOnlyKinFitLight2_[iCombi];}
  float hadrAndLeptKinFitChiSq(int iCombi) { return hadrAndLeptKinFitChiSq_[iCombi]; }
  TLorentzVector hadrAndLeptKinFitLepton(int iCombi) {return hadrAndLeptKinFitLepton_[iCombi];}
  TLorentzVector hadrAndLeptKinFitNeutrino(int iCombi) {return hadrAndLeptKinFitNeutrino_[iCombi];}
  TLorentzVector hadrAndLeptKinFitLeptB(int iCombi) {return hadrAndLeptKinFitLeptB_[iCombi];}
  TLorentzVector hadrAndLeptKinFitHadrB(int iCombi) {return hadrAndLeptKinFitHadrB_[iCombi];}
  TLorentzVector hadrAndLeptKinFitLight1(int iCombi) {return hadrAndLeptKinFitLight1_[iCombi];}
  TLorentzVector hadrAndLeptKinFitLight2(int iCombi) {return hadrAndLeptKinFitLight2_[iCombi];}
  int hadrBJet() const { return hadrBJet_; } // index, according to MC
  int hadrLJet1() const { return hadrLJet1_; } // index, according to MC
  int hadrLJet2() const { return hadrLJet2_; } // index, according to MC
  int leptBJet() const { return leptBJet_; } // index, according to MC
  TLorentzVector met() const { return MET_; }
  vector<TLorentzVector> selectedJets() const { return selectedJets_; }
  TLorentzVector lepton() const {return lepton_;}
  vector<float> bTagTCHE() const { return bTagTCHE_; }
  vector<float> bTagTCHP() const { return bTagTCHP_; }
  vector<float> bTagSSVHE() const { return bTagSSVHE_; }
  vector<float> bTagSSVHP() const { return bTagSSVHP_; }
  vector<float> bTagCSV() const {return bTagCSV_;}
  vector<float> bTagJP() const {return bTagJP_;}
  vector<float> bTagJBP() const {return bTagJBP_;}
  TLorentzVector hadrBQuark() const { return hadrBQuark_; }
  TLorentzVector hadrLQuark1() const { return hadrLQuark1_; }
  TLorentzVector hadrLQuark2() const { return hadrLQuark2_; }
  TLorentzVector leptBQuark() const { return leptBQuark_; }
  float standardCosTheta() const {return standardCosTheta_;}
  TLorentzVector standardNeutrino() const {return standardNeutrino_;}
  TLorentzVector standardLepton() const {return standardLepton_;}
  float recMassW() const { return recMassW_;}
  float recMassTop() const { return recMassTop_;}
  float recMassWLept() const { return recMassWLept_;}
  float recMassTopLept() const { return recMassTopLept_;}
  float relIso() const { return relIso_;}
  int isoMuTriggerBool() const {return isoMuTriggerBool_;}
  int muonTriggerValue() const {return muonTriggerValue_;}
  int electronTriggerValue() const {return electronTriggerValue_;}
  float deltaRJetLepton() const {return deltaRJetLepton_;}
  float deltaRMuonJet() const {return deltaRMuonJet_;}
  float deltaRElectronJet() const {return deltaRElectronJet_;}
  
  void setEventID(unsigned int eventID) { eventID_ = eventID; }
  void setRunID(unsigned int runID) { runID_ = runID; }
  void setLumiBlockID(unsigned int lumiBlockID) { lumiBlockID_ = lumiBlockID; }
  void setNPV(unsigned int nPV) { nPV_ = nPV; }
  void setNPUBXm1(unsigned int nPUBXm1) { nPUBXm1_ = nPUBXm1; }
  void setNPU(unsigned int nPU) { nPU_ = nPU; }
  void setNPUBXp1(unsigned int nPUBXp1) { nPUBXp1_ = nPUBXp1; }
  //  void setEventWeight(float eventWeight) { eventWeight_ = eventWeight; }
  void setKinFitHadr(vector<float> HadrKinFitChiSq, vector<TLorentzVector> HadrKinFitHadrB, vector<TLorentzVector> HadrKinFitLight1, vector<TLorentzVector> HadrKinFitLight2){
    for(unsigned int i=0; i<12; i++){
      hadrKinFitChiSq_[i] = HadrKinFitChiSq[i];
     hadrKinFitHadrB_[i] = HadrKinFitHadrB[i];
      hadrKinFitLight1_[i] = HadrKinFitLight1[i];
      hadrKinFitLight2_[i] = HadrKinFitLight2[i];
    }
  }
  void setKinFitHadrAndLeptWOnly(vector<float> HadrAndLeptWOnlyKinFitChiSq, vector<TLorentzVector> HadrAndLeptWOnlyKinFitLepton, vector<TLorentzVector> HadrAndLeptWOnlyKinFitNeutrino, vector<TLorentzVector> HadrAndLeptWOnlyKinFitHadrB, vector<TLorentzVector> HadrAndLeptWOnlyKinFitLight1, vector<TLorentzVector> HadrAndLeptWOnlyKinFitLight2){
    for(unsigned int i=0; i<12; i++){
      hadrAndLeptWOnlyKinFitChiSq_[i] = HadrAndLeptWOnlyKinFitChiSq[i];
      hadrAndLeptWOnlyKinFitLepton_[i] = HadrAndLeptWOnlyKinFitLepton[i];
      hadrAndLeptWOnlyKinFitNeutrino_[i] = HadrAndLeptWOnlyKinFitNeutrino[i];
      hadrAndLeptWOnlyKinFitHadrB_[i] = HadrAndLeptWOnlyKinFitHadrB[i];
      hadrAndLeptWOnlyKinFitLight1_[i] = HadrAndLeptWOnlyKinFitLight1[i];
      hadrAndLeptWOnlyKinFitLight2_[i] = HadrAndLeptWOnlyKinFitLight2[i];
    }
  }
  void setKinFitHadrAndLept(vector<float> HadrAndLeptKinFitChiSq, vector<TLorentzVector> HadrAndLeptKinFitLepton, vector<TLorentzVector> HadrAndLeptKinFitNeutrino, vector<TLorentzVector> HadrAndLeptKinFitLeptB, vector<TLorentzVector> HadrAndLeptKinFitHadrB, vector<TLorentzVector> HadrAndLeptKinFitLight1, vector<TLorentzVector> HadrAndLeptKinFitLight2){
    for(unsigned int i=0; i<12; i++){
      hadrAndLeptKinFitChiSq_[i] = HadrAndLeptKinFitChiSq[i];
      hadrAndLeptKinFitLepton_[i] = HadrAndLeptKinFitLepton[i];
      hadrAndLeptKinFitNeutrino_[i] = HadrAndLeptKinFitNeutrino[i];
      hadrAndLeptKinFitLeptB_[i] = HadrAndLeptKinFitLeptB[i];
      hadrAndLeptKinFitHadrB_[i] = HadrAndLeptKinFitHadrB[i];
      hadrAndLeptKinFitLight1_[i] = HadrAndLeptKinFitLight1[i];
      hadrAndLeptKinFitLight2_[i] = HadrAndLeptKinFitLight2[i];
    }
  }
  void setHadrBJet(int hadrBJet) { hadrBJet_ = hadrBJet; }
  void setHadrLJet1(int hadrLJet1) { hadrLJet1_ = hadrLJet1; }
  void setHadrLJet2(int hadrLJet2) { hadrLJet2_ = hadrLJet2; }
  void setLeptBJet(int leptBJet) { leptBJet_ = leptBJet; }
  void setMET(TLorentzVector MET) { MET_ = MET; }
  void setSelectedJets(vector<TLorentzVector> selectedJets) { selectedJets_ = selectedJets; }
  void setLepton(TLorentzVector lepton) { lepton_ = lepton; }
  void setBTagTCHE(vector<float> bTagTCHE) { bTagTCHE_ = bTagTCHE; }
  void setBTagTCHP(vector<float> bTagTCHP) { bTagTCHP_ = bTagTCHP; }
  void setBTagSSVHE(vector<float> bTagSSVHE) { bTagSSVHE_ = bTagSSVHE; }
  void setBTagSSVHP(vector<float> bTagSSVHP) { bTagSSVHP_ = bTagSSVHP; }
  void setBTagCSV(vector<float> bTagCSV) {bTagCSV_ = bTagCSV;}
  void setBTagJP(vector<float> bTagJP) {bTagJP_ = bTagJP;}
  void setBTagJBP(vector<float> bTagJBP) {bTagJBP_ = bTagJBP;}
  void setHadrBQuark(TLorentzVector hadrBQuark) { hadrBQuark_ = hadrBQuark; }
  void setHadrLQuark1(TLorentzVector hadrLQuark1) { hadrLQuark1_ = hadrLQuark1; }
  void setHadrLQuark2(TLorentzVector hadrLQuark2) { hadrLQuark2_ = hadrLQuark2; }
  void setLeptBQuark(TLorentzVector leptBQuark) { leptBQuark_ = leptBQuark; }
  void setStandardCosTheta(float standardCosTheta) { standardCosTheta_ = standardCosTheta;}
  void setStandardNeutrino(TLorentzVector standardNeutrino) { standardNeutrino_ = standardNeutrino;}
  void setStandardLepton(TLorentzVector standardLepton) { standardLepton_ = standardLepton;}
  void setCorrectRecMassW( float recMassW) {recMassW_ = recMassW; }
  void setCorrectRecMassTop( float recMassTop) {recMassTop_ = recMassTop;}
  void setCorrectRecMassWLept( float recMassWLept) { recMassWLept_ = recMassWLept;}
  void setCorrectRecMassTopLept( float recMassTopLept) { recMassTopLept_ = recMassTopLept;}
  void setRelativeIsolation( float relIso) { relIso_ = relIso;}
  void setIsoMuTriggerBool( int isoMuTriggerBool) { isoMuTriggerBool_ = isoMuTriggerBool;}
  void setMuonTriggerValue( int muonTriggerValue) { muonTriggerValue_ = muonTriggerValue;}
  void setElectronTriggerValue( int electronTriggerValue) { electronTriggerValue_ = electronTriggerValue;}
  void setDeltaRJetLepton( float deltaRJetLepton) { deltaRJetLepton_ = deltaRJetLepton;}
  void setDeltaRMuonJet( float deltaRMuonJet) { deltaRMuonJet_ = deltaRMuonJet;}
  void setDeltaRElectronJet( float deltaRElectronJet) { deltaRElectronJet_ = deltaRElectronJet;}
  
 protected:
  
  unsigned int eventID_;
  unsigned int runID_;
  unsigned int lumiBlockID_;
  unsigned int nPV_;
  unsigned int nPUBXm1_;
  unsigned int nPU_;
  unsigned int nPUBXp1_;
  //  float eventWeight_;
  float hadrKinFitChiSq_[12];
  TLorentzVector hadrKinFitHadrB_[12];
  TLorentzVector hadrKinFitLight1_[12];
  TLorentzVector hadrKinFitLight2_[12];
  float hadrAndLeptWOnlyKinFitChiSq_[12];
  TLorentzVector hadrAndLeptWOnlyKinFitLepton_[12];
  TLorentzVector hadrAndLeptWOnlyKinFitNeutrino_[12];
  TLorentzVector hadrAndLeptWOnlyKinFitHadrB_[12];
  TLorentzVector hadrAndLeptWOnlyKinFitLight1_[12];
  TLorentzVector hadrAndLeptWOnlyKinFitLight2_[12];
  float hadrAndLeptKinFitChiSq_[12];
  TLorentzVector hadrAndLeptKinFitLepton_[12];
  TLorentzVector hadrAndLeptKinFitNeutrino_[12];
  TLorentzVector hadrAndLeptKinFitLeptB_[12];
  TLorentzVector hadrAndLeptKinFitHadrB_[12];
  TLorentzVector hadrAndLeptKinFitLight1_[12];
  TLorentzVector hadrAndLeptKinFitLight2_[12];
  int hadrBJet_; //index according to MC
  int hadrLJet1_;
  int hadrLJet2_;
  int leptBJet_;
  TLorentzVector MET_;
  vector<TLorentzVector> selectedJets_; // all selected jets
  TLorentzVector lepton_;
  vector<float> bTagTCHE_; // indices like selectedJets indices
  vector<float> bTagTCHP_;
  vector<float> bTagSSVHE_;
  vector<float> bTagSSVHP_;
  vector<float> bTagCSV_;
  vector<float> bTagJP_;
  vector<float> bTagJBP_;
  TLorentzVector hadrBQuark_;
  TLorentzVector hadrLQuark1_;
  TLorentzVector hadrLQuark2_;
  TLorentzVector leptBQuark_;
  float standardCosTheta_;
  TLorentzVector standardNeutrino_;
  TLorentzVector standardLepton_;
  float recMassW_;
  float recMassTop_;
  float recMassWLept_;
  float recMassTopLept_;
  float relIso_;
  int isoMuTriggerBool_;
  int muonTriggerValue_;
  int electronTriggerValue_;
  float deltaRJetLepton_;
  float deltaRMuonJet_;
  float deltaRElectronJet_;
  
  ClassDef (WTree,2);
};

#endif

