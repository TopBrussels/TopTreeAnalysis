#ifndef VLQTree_h
#define VLQTree_h

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
#include "TLorentzVector.h"

using namespace std;

class VLQTree : public TObject
{
 public:
	
  VLQTree() :
    TObject()
    ,eventID_(0)
    ,runID_(0)
    ,lumiBlockID_(0)
    ,nPV_(0)
    ,nTruePU_(0)
		,MuTrigged_(0)
		,ElTrigged_(0)
    ,SelectednLeptons_(0)
		,SelectednMu_(0)
		,SelectednEl_(0)
		,eventWeight_(0)
		,MET_()
    ,selectedJets_()
    ,selectedJets_bTagCSV_()
		,selectedJets_partonFlavour_()
		,selectedForwardJets_()
		,selectedForwardJets_partonFlavour_()
    ,selectedMuons_()
		,selectedElectrons_()
		,selectedMuonsRelIso_()
		,selectedElectronsRelIso_()
		,selectedMuonsCharge_(9999)
		,selectedElectronsCharge_(9999)
		//,chargeMisIdRateBarrel_(0)
		//,chargeMisIdRateEndcap_(0)
		//,mcQuarksForMatching_()
   {;}
  
  ~VLQTree() {;}

  unsigned int eventID() const { return eventID_; }
  unsigned int runID() const { return runID_; }
  unsigned int lumiBlockID() const { return lumiBlockID_; }
  unsigned int nPV() const { return nPV_; }
  unsigned int nTruePU() const { return nTruePU_; }
	bool MuTrigged() const { return MuTrigged_; }
	bool ElTrigged() const { return ElTrigged_; }
  unsigned int SelectednLeptons() const { return SelectednLeptons_; }
	unsigned int SelectednMu() const { return SelectednMu_; }
	unsigned int SelectednEl() const { return SelectednEl_; }	
  float eventWeight() { return eventWeight_; }
	TLorentzVector met() const { return MET_; }  //to be adapted, only float!!
  vector<TLorentzVector> selectedJets() const { return selectedJets_; }
  vector<float> selectedJets_bTagCSV() const { return selectedJets_bTagCSV_; }
	vector<int> selectedJets_partonFlavour() const { return selectedJets_partonFlavour_; }
	vector<TLorentzVector> selectedForwardJets() const { return selectedForwardJets_; }
	vector<int> selectedForwardJets_partonFlavour() const { return selectedForwardJets_partonFlavour_; }
	vector<TLorentzVector> selectedMuons() const { return selectedMuons_; }
  vector<TLorentzVector> selectedElectrons() const { return selectedElectrons_; }
	vector<float> selectedMuonsRelIso() const { return selectedMuonsRelIso_; }
  vector<float> selectedElectronsRelIso() const { return selectedElectronsRelIso_; }
	vector<int> selectedMuonsCharge() const { return selectedMuonsCharge_; }
  vector<int> selectedElectronsCharge() const { return selectedElectronsCharge_; }
  //float chargeMisIdRateBarrel() { return chargeMisIdRateBarrel_; }
  //float chargeMisIdRateEndcap() { return chargeMisIdRateEndcap_; }
  //vector<TLorentzVector> mcQuarksForMatching() const { return mcQuarksForMatching_; }

	void setEventID(unsigned int eventID) { eventID_ = eventID; }
  void setRunID(unsigned int runID) { runID_ = runID; }
  void setLumiBlockID(unsigned int lumiBlockID) { lumiBlockID_ = lumiBlockID; }  
  void setNPV(unsigned int nPV) { nPV_ = nPV; }
  void setNTruePU(unsigned int nTruePU) { nTruePU_ = nTruePU; }
	void setMuTrigged(bool MuTrigged) { MuTrigged_ = MuTrigged; }
	void setElTrigged(bool ElTrigged) { ElTrigged_ = ElTrigged; }
	void setSelectednLeptons(unsigned int SelectednLeptons) { SelectednLeptons_ = SelectednLeptons; }
	void setSelectednMu(unsigned int SelectednMu) { SelectednMu_ = SelectednMu; }
	void setSelectednEl(unsigned int SelectednEl) { SelectednEl_ = SelectednEl; }
	void setEventWeight(float eventWeight) { eventWeight_ = eventWeight; }
	void setMET(TLorentzVector MET) { MET_ = MET; }
	void setSelectedJets(vector<TLorentzVector> selectedJets) { selectedJets_ = selectedJets; }
  void setSelectedJets_bTagCSV(vector<float> selectedJets_bTagCSV) { selectedJets_bTagCSV_ = selectedJets_bTagCSV; }
	void setSelectedJets_partonFlavour(vector<int> selectedJets_partonFlavour) { selectedJets_partonFlavour_ = selectedJets_partonFlavour; }
	void setSelectedForwardJets(vector<TLorentzVector> selectedForwardJets) { selectedForwardJets_ = selectedForwardJets; }
	void setSelectedForwardJets_partonFlavour(vector<int> selectedForwardJets_partonFlavour) { selectedForwardJets_partonFlavour_ = selectedForwardJets_partonFlavour; }
	void setMuons(vector<TLorentzVector> selectedMuons) { selectedMuons_ = selectedMuons; }
	void setElectrons(vector<TLorentzVector> selectedElectrons) { selectedElectrons_ = selectedElectrons; }
	void setMuonsRelIso(vector<float> selectedMuonsRelIso) { selectedMuonsRelIso_ = selectedMuonsRelIso; }
	void setElectronsRelIso(vector<float> selectedElectronsRelIso) { selectedElectronsRelIso_ = selectedElectronsRelIso; }
	void setMuonsCharge(vector<int> selectedMuonsCharge) { selectedMuonsCharge_ = selectedMuonsCharge; }
	void setElectronsCharge(vector<int> selectedElectronsCharge) { selectedElectronsCharge_ = selectedElectronsCharge; }	
	//void setChargeMisIdRateBarrel(float chargeMisIdRateBarrel) { chargeMisIdRateBarrel_ = chargeMisIdRateBarrel; }	
	//void setChargeMisIdRateEndcap(float chargeMisIdRateEndcap) { chargeMisIdRateEndcap_ = chargeMisIdRateEndcap; }	
	//void setmcQuarksForMatching(vector<TLorentzVector> mcQuarksForMatching) { mcQuarksForMatching_ = mcQuarksForMatching; }

 protected:
  
  unsigned int eventID_;
  unsigned int runID_;
  unsigned int lumiBlockID_;
  unsigned int nPV_;
  unsigned int nTruePU_;
	bool MuTrigged_;
	bool ElTrigged_;
  unsigned int SelectednLeptons_;
  unsigned int SelectednMu_;
	unsigned int SelectednEl_;
	float eventWeight_;
	TLorentzVector MET_;
	vector<TLorentzVector> selectedJets_; // all selected jet
	vector<float> selectedJets_bTagCSV_; // indices like selectedJets indices
	vector<int> selectedJets_partonFlavour_;
	vector<TLorentzVector> selectedForwardJets_; // all selected forward jet
	vector<int> selectedForwardJets_partonFlavour_;
	vector<TLorentzVector> selectedMuons_;
	vector<TLorentzVector> selectedElectrons_;
	vector<float> selectedMuonsRelIso_;
	vector<float> selectedElectronsRelIso_;	
	vector<int> selectedMuonsCharge_;
	vector<int> selectedElectronsCharge_;	
  //float chargeMisIdRateBarrel_;
  //float chargeMisIdRateEndcap_;
	//vector<TLorentzVector> mcQuarksForMatching_;
	
	ClassDef (VLQTree,3);
};

#endif
