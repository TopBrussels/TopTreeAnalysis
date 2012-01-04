#ifndef Monster_h
#define Monster_h

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

class Monster : public TObject
{
 public:
	
  Monster() :
    TObject()
    ,eventID_(0)
    ,runID_(0)
    ,lumiBlockID_(0)
    ,selectedSemiMu_(false)
    ,all4JetsMCMatched_(false)
    ,allHadronicJetsMCMatched_(false)
    ,mvaVals_()
    ,mvaResults_()
    ,eventWeight_(0)
    ,fitResults_()
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
    ,lepton_()
    ,leptonCharge_(0)
    ,hadrBQuark_()
    ,hadrLQuark1_()
    ,hadrLQuark2_()
    ,leptBQuark_()
   {;}
  
  ~Monster() {;}
  
  unsigned int eventID() const { return eventID_; }
  unsigned int runID() const { return runID_; }
  unsigned int lumiBlockID() const { return lumiBlockID_; }
  bool selectedSemiMu() const { return selectedSemiMu_; }
  bool all4JetsMCMatched() const { return all4JetsMCMatched_; }
  bool allHadronicJetsMCMatched() const { return allHadronicJetsMCMatched_; }
  float* mvaVals() { return mvaVals_; }
  float mvaVal(int i) const { return mvaVals_[i]; }
//  unsigned int** mvaResults() { return mvaResults_; }
  unsigned int* mvaResult(int i) { return mvaResults_[i]; }
  float eventWeight() const { return eventWeight_; }
//  vector<TH2F> fitResults() const { return TH2F(); /*fitResults_;*/ }
//  TH2F fitResult(int i) const { return TH2F(); /*fitResults_[i];*/ }
//  float*** fitResults() { return fitResults_; }
  float fitResult(int iCombi, int iJESl, int iJESb) { return fitResults_[iCombi][iJESl][iJESb]; }
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
  TLorentzVector lepton() const { return lepton_; }
  int leptonCharge() const { return leptonCharge_; }
  TLorentzVector hadrBQuark() const { return hadrBQuark_; }
  TLorentzVector hadrLQuark1() const { return hadrLQuark1_; }
  TLorentzVector hadrLQuark2() const { return hadrLQuark2_; }
  TLorentzVector leptBQuark() const { return leptBQuark_; }
  
  void setEventID(unsigned int eventID) { eventID_ = eventID; }
  void setRunID(unsigned int runID) { runID_ = runID; }
  void setLumiBlockID(unsigned int lumiBlockID) { lumiBlockID_ = lumiBlockID; }
  void setSelectedSemiMu(bool selectedSemiMu) { selectedSemiMu_ = selectedSemiMu; }
  void setAll4JetsMCMatched(bool all4JetsMCMatched) { all4JetsMCMatched_ = all4JetsMCMatched; }
  void setAllHadronicJetsMCMatched(bool allHadronicJetsMCMatched) { allHadronicJetsMCMatched_ = allHadronicJetsMCMatched; }
  void setMvaVals(vector<float> mvaVals)
  {
    for(unsigned int i=0; i<12; i++)
      mvaVals_[i] = mvaVals[i];
  }
  void setMvaResults(vector< vector<unsigned int> > mvaResults)
  {
    for(unsigned int iCombi=0; iCombi<12; iCombi++)
      for(unsigned int iJet=0; iJet<12; iJet++)
        mvaResults_[iCombi][iJet] = mvaResults[iCombi][iJet];
  }
  void setEventWeight(float eventWeight) { eventWeight_ = eventWeight; }
  void setFitResults(vector<TH2F> fitResults, bool topMass)
  {
    for(unsigned int iCombi=0; iCombi<12; iCombi++)
    {
      for(unsigned int iJESl=0; iJESl<41; iJESl++)
      {
        for(unsigned int iJESb=0; iJESb<41; iJESb++)
        {
          if( !topMass || iJESb < 26 )
            fitResults_[iCombi][iJESl][iJESb] = fitResults[iCombi].GetBinContent(iJESl+1,iJESb+1);
          else
            fitResults_[iCombi][iJESl][iJESb] = 0;
        }
      }
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
  void setLepton(TLorentzVector lepton) { lepton_ = lepton; }
  void setLeptonCharge(int leptonCharge) { leptonCharge_ = leptonCharge; }
  void setHadrBQuark(TLorentzVector hadrBQuark) { hadrBQuark_ = hadrBQuark; }
  void setHadrLQuark1(TLorentzVector hadrLQuark1) { hadrLQuark1_ = hadrLQuark1; }
  void setHadrLQuark2(TLorentzVector hadrLQuark2) { hadrLQuark2_ = hadrLQuark2; }
  void setLeptBQuark(TLorentzVector leptBQuark) { leptBQuark_ = leptBQuark; }
  
 protected:
  
  unsigned int eventID_;
  unsigned int runID_;
  unsigned int lumiBlockID_;
  bool selectedSemiMu_;
  bool all4JetsMCMatched_;
  bool allHadronicJetsMCMatched_;
  float mvaVals_[12];
  unsigned int mvaResults_[12][4]; // jet indices
  float eventWeight_;
//  vector<TH2F> fitResults_;
  float fitResults_[12][41][41];
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
  TLorentzVector lepton_;
  int leptonCharge_;
  TLorentzVector hadrBQuark_;
  TLorentzVector hadrLQuark1_;
  TLorentzVector hadrLQuark2_;
  TLorentzVector leptBQuark_;
  
  ClassDef (Monster,2);
};

#endif

